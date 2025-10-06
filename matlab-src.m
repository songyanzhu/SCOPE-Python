// src/+lut/change_detection.m
function change_detection(rmse_im_path, gpr_im_path)

    if nargin == 0
        in_dir = '../exercise';
        rmse_im_path = fullfile(in_dir, 'results_rmse.tif');
        gpr_im_path = fullfile(in_dir, 'results_gpr.tif');
    end
    out_path = fullfile(fileparts(rmse_im_path), 'results_rmse_gpr.png');
    
    assert(exist(rmse_im_path, 'file') ~= 0, ['Did not find `%s` image.\n'... 
        'Have you retrieved on full_set.csv with lut.use_rmse()?'], rmse_im_path)
    assert(exist(gpr_im_path, 'file') ~= 0, ['Did not find `%s` image.\n'... 
        'Have you retrieved on full_set.csv with lut.use_gpr()?'], gpr_im_path)
    
    lut_im = imread(rmse_im_path);
    gpr_im = imread(gpr_im_path);
    
    im_dif = lut_im - gpr_im;
    lut.plot_image(im_dif, '(RMSE - GPR)', out_path)
end

// src/+lut/compress_geotiff.m

function compress_geotiff(path_in, path_out)

    ver_out = ver;
    toolboxes = {ver_out.Name};
    geo_tif = any(strcmp('Mapping Toolbox', toolboxes));
    
    if geo_tif
        fprintf('copying georeference tags from input image %s\n', path_in)
        geoinfo = geotiffinfo(path_in);
        key = geoinfo.GeoTIFFTags.GeoKeyDirectoryTag;
        R = geoinfo.SpatialRef;
    else
        warning(['Mapping Toolbox is not installed. Output .tifs can not be georeferenced.\n'...
            'Use gdal_translate to georeference from input image %s\n'...
            'Output .tifs are identical to input image'], path_in)
    end
    
    im = imread(path_in);
    comp.Compression = Tiff.Compression.PackBits;
    geotiffwrite(path_out, im, R, 'GeoKeyDirectoryTag', key, 'TiffTags', comp)
    
end

// src/+lut/csv2image_plane.m
function im = csv2image_plane(val_in, res)
% to read r, c from column name ind_r_c
    vars = val_in.Properties.VariableNames;
    i_ind = ~cellfun(@isempty, strfind(vars, 'ind'));  % ignore readability, or check > 2016b
    splt =  strsplit(vars{i_ind}, '_');
    r = str2num(splt{2});
    c = str2num(splt{3});
    im = nan(r, c);
    pix_ind = table2array(val_in(:, i_ind));
    im(pix_ind) = res;
end

// src/+lut/generate_lut_input.m

function generate_lut_input(tab, n_spectra, outdir)
    % params = generate_lut_input(tab, n_spectra, outdir)
    if nargin == 0
        tab = readtable('+lut/input_borders.csv');
        n_spectra = 1000;
        outdir = '../exercise';
    end
    out_file = fullfile(outdir, 'lut_in.csv');
    assert(exist(out_file, 'file') == 0, '`%s` file already exists, delete it first', out_file)
    
    include = logical(tab.include);
    lb = tab.lower(include)';
    ub = tab.upper(include)';
    varnames = tab.variable(include)';

    % one row - one set of parameters
    lh = lhsdesign(n_spectra, sum(include));

    params = (ub-lb) .* lh + lb;

    if any(strcmp('LIDFa' , varnames))
        % abs(LIDFa + LIDFb) <= 1
        i_lidfa = strcmp('LIDFa', varnames);
        i_lidfb = strcmp('LIDFb', varnames);
        lidfa = params(:, i_lidfa);
        lidfb = params(:, i_lidfb);
        params(:, i_lidfa) = (lidfa + lidfb) / 2;
        params(:, i_lidfb) = (lidfa - lidfb) / 2;
    end

    t = array2table(params);
    t.Properties.VariableNames = varnames;
    t.t = datestr(datetime(2022, 7, 1, 0, 12, 1:n_spectra), 'yyyymmddHHMMSS.FFF');
    writetable(t, out_file)
    
    if verLessThan('matlab', '9.1')  % < 2016b
        varnames_in = '';
    else
        varnames_in = strjoin(varnames, ', ');
    end
    fprintf('Sampled %i parameters: %s\n', length(varnames), varnames_in)
    fprintf('Saved lut input (parameters) in `%s`\n', out_file)
end

// src/+lut/image2csv.m

function image2csv(im_path, outdir)
    
    if nargin == 0
        im_path = '../exercise/images';
        outdir = '../exercise';
        name = 'full_set';
    else
        [~, name, ~] = fileparts(im_path);
    end
    csv_out = fullfile(outdir, sprintf('%s.csv', name));
    var_names = {'Cab', 'LAI'};  % or dir(images) if only .tif
    assert(exist(im_path, 'dir') == 7, '`%s` does not exist. Please, put images (Cab.tif and LAI.tif) into that folder', im_path)
    fprintf('reading images from `%s`\n', im_path)
    
    [~, ~, ext] = fileparts(im_path);
    if strcmp(ext, '.nc')
        get_val = @(x) ncread(im_path, x);
    else
        get_val = @(x) imread(fullfile(im_path, sprintf('%s.tif', x)));
    end

    [r, c] = size(get_val(var_names{1}));
    
    %% flattening
    fprintf('Flattening %d rows, %d columns\n', r, c)
    vals = nan(r*c, length(var_names));
    for i=1:length(var_names)
        var = var_names{i};
        v = get_val(var);
%         v = v(ir, ic);
        vals(:, i) = v(:);  % flattening
    end
    
    i_nans = any(isnan(vals), 2);
    df = array2table(vals, 'VariableNames', var_names);
    df_clean = df(~i_nans, :);
    df_clean.(['ind_' num2str(r), '_', num2str(c)]) = find(~i_nans);  % or fprintf()
    
    fprintf('found %d not-nan pixels\n', sum(~i_nans))
    
    %% saving
%     out = struct();
%     out.i_nans = i_nans;
%     out.df_clean = df_clean;
%     out.r = r;
%     out.c = c;
% %     reshape(i_nans, r, c);
%     save(fullfile(outdir, sprintf('%s.mat', name)), 'out')
    
    writetable(df_clean, csv_out)
    fprintf('saved to `%s`\n', csv_out)
end

// src/+lut/input_borders.csv

variable	lower	upper	units	include	default	meaning
Cab	0	100	ug cm-2	1	80	leaf chlorophyll content
Cca	0	25	ug cm-2	1	20	leaf carotenoid content
Cant	0	5	ug cm-2	1	0	leaf anthocyanin content
Cdm	0	0.02	g cm-2	1	0.012	leaf mass per area (dry matter)
Cw	0	0.2	cm	1	0.009	equivalent leaf water thickness
Cs	0	1.2	-	1	0	senescent material (brown pigments)
N	1	3.5	-	1	1.4	mesophyll structure parameter
LAI	0	7	m2 m-2	1	3	canopy leaf area index
LIDFa	-1	1	-	1	-0.35	leaf inclination distribution function parameters a, b
LIDFb	-1	1	-	1	-0.15	
B	0	0.9	-	0	0.5	soil brightness
BSMlat	20	40	-	0	25	BSM model parameter lat
BSMlon	40	60	-	0	45	BSM model parameter lon
SMC	5	55	%	0	25	volumetric soil moisture content
Vcmax25	1	200	umol m-2 s-1	1	70	maximum carboxylation capacity (at optimum temperature)
BallBerrySlope	5	20	-	1	12	Ball-Berry stomatal conductance parameter
Rin	5	800	W m-2	1	600	broadband incoming shortwave radiation (0.4-2.5 um)
Rli	200	500	W m-2	0	300	broadband incoming longwave radiation (2.5-50 um)
Ta	5	35	C	1	20	air temperature
p	900	1100	hPa	0	970	air pressure
ea	0	60	hPa	0	15	atmospheric vapour pressure
RH	0	1	-	1	0.64	relative humidity (proxy of ea)
u	0.5	10	m s-1	0	2	wind speed at height z
tts	0	60	deg	1	30	solar zenith angle
tto	0	60	deg	1	0	observation zenith angle
psi	0	180	deg	1	90	azimuthal difference between solar and observation angle

// src/+lut/lut_search.m

function [res, res_std] = lut_search(params, lut_params, response)

    %% check validity: all lut params in params
    p_names = params.Properties.VariableNames;
    lut_names = lut_params.Properties.VariableNames;

    % absent = setdiff(lut_names, p_names);
    % assert(isempty(absent), '%s parameter must be in input table, because it is in lut\n', absent)
    % 
    % extra = setdiff(p_names, lut_names);
    % if ~isempty(extra)
    %     warning('%s parameters are not in lut, thus will not effect\n', strjoin(extra, ','))
    % end

    %% normalization to standardize RMSE, column alignment
    p_ordered = nan(size(params, 1), length(lut_names));
    for i=1:length(lut_names)
        v = lut_names{i};
        v_min = min(lut_params.(v));
        v_max = max(lut_params.(v));
        lut_params.(v) = (lut_params.(v) - v_min) / (v_max - v_min);
        p_ordered(:, i) = (params.(v) - v_min) / (v_max - v_min);
    end
    lut_arr = table2array(lut_params);

    tic
    ind = arrayfun(@(x) top_indices(x, p_ordered, lut_arr), 1:size(p_ordered, 1), 'UniformOutput', false);
    % rmses = arrayfun(@(x) top_indices(x, p_ordered, lut_arr), 1:size(p_ordered, 1), 'UniformOutput', false);
    ind = cell2mat(ind');
    toc

    actot = response(ind);

    res = nanmedian(actot, 2);
    res_std = nanstd(actot, 0, 2);  % (v, flag, dim)

end


function ind = top_indices(i, meas_arr, lut_arr)

    if mod(i, 100000) == 0
        fprintf('done %i / %i  pixels\n', i, length(meas_arr))
    end
    
    rmses = sqrt(mean((meas_arr(i, :) - lut_arr) .^ 2, 2));
    
    [~, I] = sort(rmses);
    ind = I(1:10)';  % top 10
%     rmses = rmses(ind);
end

// src/+lut/pick100.m

function pick100(full_set_path, n_pixels)
    if nargin == 0
        full_set_path = '../exercise/full_set.csv';
        n_pixels = 100;
    end
    out_file = fullfile(fileparts(full_set_path), 'validation_in.csv');
    
    assert(exist(out_file, 'file') == 0, '`%s` file already exists, delete it first', out_file)
    
    df = readtable(full_set_path);
    ind = randi(size(df, 1), n_pixels, 1);
    df_out = df(ind, :);
    
    writetable(df_out, out_file)
    fprintf('saved 100 pixels to `%s`\n', out_file)
end

// src/+lut/plot_1to1.m

function plot_1to1(meas, mod, what, out_path, var_name)
    if nargin < 5
        var_name = 'GPP';
    end

    figure
%     errorbar(val_out, res, res_std, 'o')
    units = '[W m^{-2}]';
    if any(strcmp(var_name, {'Actot', 'GPP'}))
        units = '[\mumol CO_2 m^{-2} s^{-1}]';
    end
    scatter(meas, mod)
    xlabel(sprintf('SCOPE %s %s', var_name, units))
    ylabel(sprintf('%s %s %s', what, var_name, units))
    l = refline(1, 0);
    l.Color = 'r';
    refline
    rmse = sqrt(mean((meas - mod) .^ 2));
    bias = mean(mod-meas);
    lm = fitlm(meas, mod); 
    title(sprintf('%s, RMSE = %.2f, bias=%.2f, R^2_{adj} = %.2f', var_name, rmse, bias, lm.Rsquared.Adjusted))
    legend('scatter', '1-to-1 line', 'regression', 'Location', 'southeast')
    saveas(gcf, out_path)
end

// src/+lut/plot_image.m

function plot_image(im, what, out_path)
    figure
    imagesc(im)
    cb = colorbar;
    title(cb, '[\mumol CO_2 m^{-2} s^{-1}]')
    title(['GPP with ' what])
    % if GPR is out of range (for SNAP), highly negative values appear
    caxis([quantile(im(:), 0.01), quantile(im(:), 0.99)])  
    saveas(gcf, out_path)
end

// src/+lut/train_gpr.m

function train_gpr(lut_in_path, lut_out_path, var_name)
    if nargin == 0
        in_dir = '../exercise';
        lut_in_path = fullfile(in_dir, 'lut_in.csv');
        lut_out_path = fullfile(in_dir, 'lut_out.csv');
        var_name = 'Actot';
    end
    mat_out = fullfile(fileparts(lut_in_path), ['gpr_' var_name '.mat']);
    
    assert(exist(lut_in_path, 'file') ~= 0, ['Did not find `%s` file.\n'... 
        'Have you generated the LUT input with lut.generate_lut_input()?'], lut_in_path)
    assert(exist(lut_out_path, 'file') ~= 0, ['Did not find `%s` file.\n'...
        'Have you copied the `fluxes.csv` from `../output` after SCOPE run on lut_in.csv?'], lut_out_path)
    
    lut_in = readtable(lut_in_path);
    
    opt = detectImportOptions(lut_out_path);
    flu = readtable(lut_out_path, opt);
    lut_in.(var_name) = flu.(var_name);
    
    % Cross varidation (train: 70%, test: 30%)
    cv = cvpartition(size(lut_in,1),'HoldOut',0.3);
    idx = cv.test;
    % Separate to training and test data
    lut_train = lut_in(~idx,:);
    lut_test  = lut_in(idx,:);
    
    fprintf('Started training gaussian process regression.\nUsually takes 1 minute.\n')
    gprMdl = fitrgp(lut_train, var_name, 'Standardize', 1);
    % we can overfit but not predict on new data with cvGPR (kfoldGPR)
%     gprMdl = fitrgp(lut_train, var_name, 'Standardize', 1, 'CrossVal', 'on', 'Verbose', 1);
%     res = kfoldPredict(gprMdl);
    save(mat_out, 'gprMdl')
    fprintf('GPR is saved in `%s`\n', mat_out)
    

    res = predict(gprMdl, lut_test);
    fig_path = fullfile(fileparts(lut_in_path), [var_name '_gpr.png']);
    lut.plot_1to1(lut_test.(var_name), res, 'GPR', fig_path, var_name)
    
end

// src/+lut/use_gpr.m

function use_gpr(gpr_path, full_set_path)

    if nargin == 0
        in_dir = '../exercise';
        gpr_path = fullfile(in_dir, 'gpr.mat');
        full_set_path = fullfile(in_dir, 'full_set.csv');
    end
    csv_out = fullfile(fileparts(full_set_path), 'results_gpr.csv');
    map_out = fullfile(fileparts(full_set_path), 'results_gpr.tif');
    fig_out = fullfile(fileparts(full_set_path), 'results_gpr.png');
    
    assert(exist(gpr_path, 'file') ~= 0, ['Did not find `%s` file.\n'... 
        'Have you trained the gaussian process regression (GPR) with lut.train_gpr()?'], gpr_path)
    assert(exist(full_set_path, 'file') ~= 0, ['Did not find `%s` file.\n'... 
        'Have you flattened the images with lut.image2csv()?'], full_set_path)

    gpr = load(gpr_path);
    gprMdl = gpr.gprMdl;
    val_in = readtable(full_set_path);
    
    fprintf('Working, usually takes around 1 minute\n')
    res = predict(gprMdl, val_in);
    
    tab = array2table(res, 'VariableNames', {'Actot'});
    writetable(tab, csv_out)
    fprintf('saved `%s`\n', csv_out)
    
    im = lut.csv2image_plane(val_in, res);
    lut.write_tiff(im, map_out)
    
    lut.plot_image(im, 'GPR', fig_out)
end

// src/+lut/use_rmse.m

function use_rmse(lut_in_path, lut_out_path, full_set_path, var_name)

    if nargin == 0
        in_dir = '../exercise';
        lut_in_path = fullfile(in_dir, 'lut_in.csv');
        lut_out_path = fullfile(in_dir, 'lut_out.csv');
        full_set_path = fullfile(in_dir, 'full_set.csv');
        var_name = 'Actot';
    end
    csv_out = fullfile(fileparts(full_set_path), 'results_rmse.csv');
    map_out = fullfile(fileparts(full_set_path), 'results_rmse.tif');
    fig_out = fullfile(fileparts(full_set_path), 'results_rmse.png');
    
    assert(exist(lut_in_path, 'file') ~= 0, ['Did not find `%s` file.\n'... 
        'Have you generated the LUT input with lut.generate_lut_input()?'], lut_in_path)
    assert(exist(lut_out_path, 'file') ~= 0, ['Did not find `%s` file.\n'...
        'Have you copied the `fluxes.csv` from `../output` after SCOPE run on lut_in.csv?'], lut_out_path)
    assert(exist(full_set_path, 'file') ~= 0, ['Did not find `%s` file.\n'... 
        'Have you flattened the images with lut.image2csv()?'], full_set_path)

    lut_in = readtable(lut_in_path);
    opt = detectImportOptions(lut_out_path);
    flu = readtable(lut_out_path, opt);
    lut_out = flu.(var_name);
    val_in = readtable(full_set_path);
    
    [res, res_std] = lut.lut_search(val_in, lut_in, lut_out);
    
    tab = array2table([res, res_std], 'VariableNames', {'Actot', 'Actot_sd'});
    writetable(tab, csv_out)
    fprintf('saved `%s`\n', csv_out)
    
    im = lut.csv2image_plane(val_in, res);
    lut.write_tiff(im, map_out)
    
    lut.plot_image(im, 'RMSE', fig_out)
end

// src/+lut/validate_gpr.m

function validate_gpr(gpr_path, val_in_path, val_out_path)

    if nargin == 0
        in_dir = '../exercise';
        gpr_path = fullfile(in_dir, 'gpr.mat');
        val_in_path = fullfile(in_dir, 'validation_in.csv');
        val_out_path = fullfile(in_dir, 'validation_out.csv');
    end
    fig_path = fullfile(fileparts(val_in_path), 'validation_gpr.png');
    
    assert(exist(gpr_path, 'file') ~= 0, ['Did not find `%s` file.\n'... 
        'Have you trained the gaussian process regression (GPR) with lut.train_gpr()?'], gpr_path)
    assert(exist(val_in_path, 'file') ~= 0, ['Did not find `%s` file.\n'... 
        'Have you generated the validation input with lut.pick100()?'], val_in_path)
    assert(exist(val_out_path, 'file') ~= 0, ['Did not find `%s` file.\n'...
        'Have you copied the `fluxes.csv` from `../output` after SCOPE run on validation_in.csv?'], val_out_path)

    gpr = load(gpr_path);
    gprMdl = gpr.gprMdl;
    val_in = readtable(val_in_path);
    val_out = lut.read_actot(val_out_path);
    
    tic
    res = predict(gprMdl, val_in);
    toc
    
    lut.plot_1to1(val_out, res, 'GPR', fig_path)

end

// src/+lut/validate_rmse.m

function validate_rmse(lut_in_path, lut_out_path, val_in_path, val_out_path)

    if nargin == 0
        in_dir = '../exercise';
        lut_in_path = fullfile(in_dir, 'lut_in.csv');
        lut_out_path = fullfile(in_dir, 'lut_out.csv');
        val_in_path = fullfile(in_dir, 'validation_in.csv');
        val_out_path = fullfile(in_dir, 'validation_out.csv');
    end
    fig_path = fullfile(fileparts(val_in_path), 'validation_rmse.png');
    
    assert(exist(lut_in_path, 'file') ~= 0, ['Did not find `%s` file.\n'... 
        'Have you generated the LUT input with lut.generate_lut_input()?'], lut_in_path)
    assert(exist(lut_out_path, 'file') ~= 0, ['Did not find `%s` file.\n'...
        'Have you copied the `fluxes.csv` from `../output` after SCOPE run on lut_in.csv?'], lut_out_path)
    assert(exist(val_in_path, 'file') ~= 0, ['Did not find `%s` file.\n'... 
        'Have you generated the validation input with lut.pick100()?'], val_in_path)
    assert(exist(val_out_path, 'file') ~= 0, ['Did not find `%s` file.\n'...
        'Have you copied the `fluxes.csv` from `../output` after SCOPE run on validation_in.csv?'], val_out_path)

    lut_in = readtable(lut_in_path);
    lut_out = lut.read_actot(lut_out_path);
    val_in = readtable(val_in_path);
    val_out = lut.read_actot(val_out_path);
    
    [res, res_std] = lut.lut_search(val_in, lut_in, lut_out);
    
    lut.plot_1to1(val_out, res, 'RMSE', fig_path)

end

// src/+lut/write_tiff.m

function write_tiff(im, out_path)
    
    input_image_path = '../exercise/images/Cab.tif';
    copyfile(input_image_path, out_path)
    t = Tiff(out_path, 'r+');
    if getTag(t, 'BitsPerSample') == 32
        im = single(im);
    end
    write(t, im)
    close(t)

end

// src/IO/assignvarnames.m

function V = assignvarnames()
V                           = struct('Name','','Val', zeros(64,1));
V(1).Name                   = 'Cab';
V(2).Name                   = 'Cca';
V(3).Name                   = 'Cdm';
V(4).Name                   = 'Cw';
V(5).Name                   = 'Cs';
V(6).Name                   = 'N';
V(7).Name                  = 'Cant';       %Added March 2017
V(8).Name                   = 'Cp';
V(9).Name                  = 'Cbc'; % carbon based consituents PROSPECT-PRO
V(10).Name                   = 'rho_thermal';
V(11).Name                   = 'tau_thermal';

V(12).Name                   = 'Vcmax25';
V(13).Name                  = 'BallBerrySlope';  % see # 64, below for intercept: 'BallBerry0'
V(14).Name                  = 'BallBerry0';
V(15).Name                  = 'Type';
V(16).Name                  = 'kV';
V(17).Name                  = 'Rdparam';
V(18).Name                  = 'fqe';
V(19).Name                  = 'Kn0';
V(20).Name                  = 'Knalpha';
V(21).Name                  = 'Knbeta';

V(22).Name                  = 'LAI';
V(23).Name                  = 'hc';
V(24).Name                  = 'zo';
V(25).Name                  = 'd';
V(26).Name                  = 'LIDFa';
V(27).Name                  = 'LIDFb';
V(28).Name                  = 'leafwidth';

V(29).Name                  = 'z';
V(30).Name                  = 'Rin';
V(31).Name                  = 'Ta';
V(32).Name                  = 'Rli';
V(33).Name                  = 'p';
V(34).Name                  = 'ea';
V(35).Name                  = 'u';
V(36).Name                  = 'Ca';
V(37).Name                  = 'Oa';

V(38).Name                  = 'rb';
V(39).Name                  = 'Cd';
V(40).Name                  = 'CR';
V(41).Name                  = 'CD1';
V(42).Name                  = 'Psicor';
V(43).Name                  = 'CSSOIL';
V(44).Name                  = 'rbs';
V(45).Name                  = 'rwc';

V(46).Name                  = 'startDOY';
V(47).Name                  = 'endDOY';
V(48).Name                  = 'LAT';
V(49).Name                  = 'LON';
V(50).Name                  = 'timezn';
V(51).Name                  = 'tts';
V(52).Name                  = 'tto';
V(53).Name                  = 'psi';

V(54).Name                  = 'SMC';
V(55).Name                  = 'Tyear';
V(56).Name                  = 'beta';
V(57).Name                  = 'kNPQs';
V(58).Name                  = 'qLs';
V(59).Name                  = 'stressfactor';

V(60).Name                  = 'spectrum';
V(61).Name                  = 'BSMBrightness';
V(62).Name                  = 'BSMlat';
V(63).Name                  = 'BSMlon';

V(64).Name                  = 'rss';
V(65).Name                  = 'rs_thermal';
V(66).Name                  = 'cs';
V(67).Name                  = 'rhos';
V(68).Name                  = 'lambdas';

%V(69).Name                  = 'Cv'; 
%V(70).Name                  = 'crowndiameter'; 


end

// src/IO/bin_to_csv.m

function bin_to_csv(fnames, V, vmax, n_col, ns)

%% pars
if sum(vmax>1)
  write_output(['n_pars', {V(vmax>1).Name}], {''}, fnames.pars_file, n_col.pars, ns)
end

%% aPAR
apar_names = {'simulation_number', 'year', 'DoY', 'iPAR', 'iPARE', 'LAIsunlit', 'LAIshaded'...
    'aPARtot', 'aPARsun', 'aPARsha',...
    'aPARCabtot', 'aPARCabsun', 'aPARCabsha',...
    'aPARCartot', 'aPARCarsun', 'aPARCarsha',...
    'aPARtotE', 'aPARsunE', 'aPARshaE',...
    'aPARCabtotE', 'aPARCabsunE', 'aPARCabshaE',...
    'aPARCartotE', 'aPARCarsunE', 'aPARCarshaE',};

apar_units = {'', '', '', 'umol m-2 s-1','W m-2', 'm2m-2','m2m-2',...
    'umol m-2 s-1', 'umol m-2 s-1', 'umol m-2 s-1',...
    'umol m-2 s-1', 'umol m-2 s-1', 'umol m-2 s-1',...
    'umol m-2 s-1', 'umol m-2 s-1', 'umol m-2 s-1',...
    'W m-2','W m-2','W m-2',...
    'W m-2','W m-2','W m-2',...
    'W m-2','W m-2','W m-2'};

write_output(apar_names, apar_units, fnames.apar_file, n_col.apar, ns)

%% veg
veg_names = {'simulation_number', 'year', 'DoY', 'Photosynthesis', 'Electron_transport', 'NPQ_energy',  'NPQ_photon', 'canopy_level_FQE','LST','emis', 'GPP'};
veg_units = {'', '', '', 'umol CO2 m-2 s-1', 'umol m-2 s-1', 'W m-2', 'umol m-2 s-1', 'umol photons (umol photons)-1', 'K','', 'umol CO2 m-2 s-1'};
write_output(veg_names, veg_units, fnames.veg_file, n_col.veg, ns)

%% flu
flu_names = {'simulation_number','nu_iterations', 'year','DoY',...
    'Rnctot','lEctot','Hctot','Actot','Tcave', ...
    'Rnstot','lEstot','Hstot','Gtot','Tsave',...
    'Rntot','lEtot','Htot'};
flu_units = {'', '', '', '',  ...
    'W m-2','W m-2','W m-2','umol m-2 s-1','C',...
    'W m-2','W m-2','W m-2','W m-2','C',...
    'W m-2',' W m-2','W m-2'};
write_output(flu_names, flu_units, fnames.flu_file, n_col.flu, ns)

%% rad
rad_names = {'simulation_number','year','DoY','ShortIn','LongIn','HemisOutShort','HemisOutLong','Lo','Lot','Lote'};
rad_units = {'','','','W m-2','W m-2','W m-2','W m-2','W m-2 sr-1','W m-2 sr-1','W m-2 sr-1'};
write_output(rad_names, rad_units, fnames.rad_file, n_col.rad, ns)

%% fluor
if isfield(fnames, 'fluor_file')
    fluor_names = {'F_1stpeak', 'wl_1stpeak', 'F_2ndpeak', 'wl_2ndpeak', 'F684', 'F761', 'LFtot', 'EFtot', 'EFtot_RC'};
    fluor_units = {'W m-2 um-1 sr-1','nm','W m-2 um-1 sr-1','nm','W m-2 um-1 sr-1','W m-2 um-1 sr-1','W m-2 sr-1','W m-2','W m-2'};
    write_output(fluor_names, fluor_units, fnames.fluor_file, n_col.fluor, ns)

    write_output({'fluorescence_spectrum 640:1:850 nm'}, {'W m-2 um-1 sr-1'}, ...
        fnames.fluor_spectrum_file, n_col.fluor_spectrum, ns, true)

    write_output({'escape probability 640:1:850 nm'}, {''}, ...
        fnames.sigmaF_file, n_col.sigmaF, ns, true)

    write_output({'fluorescence_spectrum 640:1:850 nm hemispherically integrated'}, {'W m-2 um-1'}, ...
        fnames.fhemis_file, n_col.fhemis, ns, true)

    write_output({'fluorescence_spectrum 640:1:850 nm reabsorption corrected'}, {'W m-2 um-1'}, ...
        fnames.fRC_file, n_col.fhemis, ns, true)

    write_output({'fluorescence_spectrum 640:1:850 nm emission by all leaves'}, {'W m-2 um-1'}, ...
        fnames.fRCL_file, n_col.fhemis, ns, true)

    write_output({'upwelling radiance including fluorescence'}, {'W m-2 um-1 sr-1'}, ...
        fnames.Lo2_file, n_col.Lo2, ns, true)

    write_output({'apparent reflectance'}, {''}, ...
        fnames.rapp_file, n_col.rapp, ns, true)

end

%% reflectance
if isfield(fnames, 'r_file')
    write_output({'reflectance'}, {'pi*upwelling radiance/irradiance'}, ...
        fnames.r_file, n_col.r, ns, true)

    write_output({'rsd'}, {'directional-hemispherical reflectance factor'}, ...
        fnames.rsd_file, n_col.rsd, ns, true)

    write_output({'rdd'}, {'bi-hemispherical reflectance factor'}, ...
        fnames.rdd_file, n_col.rdd, ns, true)

    write_output({'rso'}, {'bi-directional reflectance factor'}, ...
        fnames.rso_file, n_col.rso, ns, true)

    write_output({'rdo'}, {'hemispherical-directional reflectance factor'}, ...
        fnames.rdo_file, n_col.rdo, ns, true)

    %% radiance
    write_output({'hemispherically integrated upwelling radiance'}, {'W m-2 um-1'}, ...
        fnames.Eout_file, n_col.Eout, ns, true)

    write_output({'upwelling radiance excluding fluorescence'}, {'W m-2 um-1 sr-1'}, ...
        fnames.Lo_file, n_col.Lo, ns, true)

    write_output({'direct solar irradiance'}, {'W m-2 um-1'}, ...
        fnames.Esun_file, n_col.Esun, ns, true)

    write_output({'diffuse solar irradiance'}, {'W m-2 um-1'}, ...
        fnames.Esky_file, n_col.Esky, ns, true)
end

%% resistances
write_output({'aerodynamicresistance(ra)', 'raforsoil', 'rss','ustar'}, ...
    {'s m-1','s m-1','s m-1','m s-1'}, ...
    fnames.resist_file, n_col.resist, ns, true)

fclose('all');

%% deleting .bin
structfun(@delete, fnames)
end

function write_output(header, units, bin_path, f_n_col, ns, not_header)
    if nargin == 5
        not_header = false;
    end
    n_csv = strrep(bin_path, '.bin', '.csv');

    f_csv = fopen(n_csv, 'w');
    header_str = [strjoin(header, ','), '\n'];
    if not_header
        header_str = ['#' header_str];
    else
        % it is a header => each column must have one
        assert(length(header) == f_n_col, 'Less headers than lines `%s` or n_col is wrong', bin_path)
    end
    fprintf(f_csv, header_str);
    fprintf(f_csv, ['#' strjoin(units, ','), '\n']);

    f_bin = fopen(bin_path, 'r');
    out = fread(f_bin, 'double');
%     fclose(f_bin);  % + some useconds to execution
    out_2d = reshape(out, f_n_col, ns)';
%     dlmwrite(n_csv, out_2d, '-append', 'precision', '%d'); % SLOW!
    for k=1:ns
        fprintf(f_csv, '%d,', out_2d(k, 1:end-1));
        fprintf(f_csv, '%d\n', out_2d(k, end));  % saves from extra comma
    end
%     fclose(f_csv);
end

// src/IO/create_output_files_binary.m

function [Output_dir, f, fnames] = create_output_files_binary(parameter_file, F, path_of_code,input_path,spectral,options)
%% Create Output dir
string          = clock;
simulation_name = char(F(1).FileName);
outdir_name     = sprintf('%s_%4.0f-%02.0f-%02.0f-%02.0f%02.0f', simulation_name, string(1:5));
Output_dir = [fullfile('output', outdir_name) filesep];

warning('off','MATLAB:DELETE:FileNotFound')
if any(~exist(Output_dir,'dir'))
    mkdir(Output_dir)
    mkdir([Output_dir,'Parameters' filesep])
end

%% Log File
for i = 1:length(parameter_file{1})
    copy_name = [strrep(parameter_file{1}{i}, '.csv', '') '_' outdir_name '.csv'];
    copyfile([input_path parameter_file{1}{i}],[Output_dir,'Parameters/', copy_name],'f')
end
fidpath          = fopen([Output_dir,'Parameters/SCOPEversion.txt'],'w');      % complete path of the SCOPE code
fprintf(fidpath,'%s', path_of_code);

%% Filenames, will become .csv if options is on

fnames.pars_file               = fullfile(Output_dir,'pars_and_input_short.bin');
fnames.apar_file               = fullfile(Output_dir,'aPAR.bin');
fnames.veg_file                = fullfile(Output_dir,'vegetation.bin');
fnames.flu_file                = fullfile(Output_dir,'fluxes.bin');
fnames.rad_file                = fullfile(Output_dir,'radiation.bin');
if options.calc_fluor
    fnames.fluor_file              = fullfile(Output_dir,'fluorescence_scalars.bin');
    fnames.fluor_spectrum_file     = fullfile(Output_dir,'fluorescence.bin');
    fnames.sigmaF_file             = fullfile(Output_dir,'sigmaF.bin');
    fnames.fhemis_file             = fullfile(Output_dir,'fluorescence_hemis.bin');
    fnames.fRC_file                = fullfile(Output_dir,'fluorescence_ReabsCorr.bin');
    fnames.fRCL_file               = fullfile(Output_dir,'fluorescence_AllLeaves.bin'); 
    fnames.Lo2_file                = fullfile(Output_dir,'Lo_spectrum_inclF.bin');
    fnames.rapp_file               = fullfile(Output_dir,'apparent_reflectance.bin');
end
if options.save_spectral
    fnames.r_file                  = fullfile(Output_dir,'reflectance.bin');
    fnames.rsd_file                = fullfile(Output_dir,'rsd.bin');
    fnames.rdd_file                = fullfile(Output_dir,'rdd.bin');
    fnames.rso_file                = fullfile(Output_dir,'rso.bin');
    fnames.rdo_file                = fullfile(Output_dir,'rdo.bin');
    fnames.Eout_file               = fullfile(Output_dir,'Eout_spectrum.bin');
    fnames.Lo_file                 = fullfile(Output_dir,'Lo_spectrum.bin');
    fnames.Esun_file               = fullfile(Output_dir,'Esun.bin');
    fnames.Esky_file               = fullfile(Output_dir,'Esky.bin');
end
fnames.resist_file             = fullfile(Output_dir,'resistances.bin');

%% Open files for writing
f = structfun(@(x) fopen(x, 'w'), fnames, 'UniformOutput',false);

%% write wl
wlS = spectral.wlS; %#ok<*NASGU>
wlF = spectral.wlF;

save([Output_dir 'wlS.txt'], 'wlS', '-ascii');
save([Output_dir 'wlF.txt'], 'wlF', '-ascii');
end

// src/IO/define_bands.m

function [spectral] = define_bands()

% Define spectral regions for SCOPE v_1.40
% All spectral regions are defined here as row vectors
% WV Jan. 2013

% 3 spectral regions for SCOPE

reg1 =   400 :    1 :  2400;
reg2 =  2500 :  100 : 15000;
reg3 = 16000 : 1000 : 50000;

spectral.wlS  = [reg1 reg2 reg3]';

% Other spectral (sub)regions  
spectral.wlP   = reg1';                            % PROSPECT data range
spectral.wlE   = 400:1:750;                       % excitation in E-F matrix
spectral.wlF   = 640:1:850;                       % chlorophyll fluorescence in E-F matrix
spectral.wlO   = reg1;                            % optical part
spectral.wlT   = [reg2 reg3];                     % thermal part
spectral.wlZ   = 500:1:600;                       % xanthophyll region
wlS            = spectral.wlS;
spectral.wlPAR = wlS(wlS>=400 & wlS<=700);  % PAR range
spectral.IwlF  = (640:850)-399;
spectral.IwlP  = 1:length(reg1);
spectral.IwlT  = length(reg1)+1: length(reg1)+length(reg2)+length(reg3);
% Data used by aggreg routine to read in MODTRAN data

spectral.SCOPEspec.nreg = 3;
spectral.SCOPEspec.start = [ 400  2500  16000];
spectral.SCOPEspec.end   = [2400 15000  50000];
spectral.SCOPEspec.res   = [   1   100   1000];

// src/IO/define_constants.m

function [const]=define_constants()

    const.A         = 6.02214E23; % [mol-1]       Constant of Avogadro
    const.h         = 6.6262E-34; % [J s]         Planck's constant
    const.c         = 299792458;  % [m s-1]       Speed of light
    const.cp        = 1004;       % [J kg-1 K-1]  Specific heat of dry air
    const.R         = 8.314;      % [J mol-1K-1]  Molar gas constant
    const.rhoa      = 1.2047;     % [kg m-3]      Specific mass of air
    const.g         = 9.81;       % [m s-2]       Gravity acceleration
    const.kappa     = 0.4;        % []            Von Karman constant
    const.MH2O      = 18;         % [g mol-1]     Molecular mass of water
    const.Mair      = 28.96;      % [g mol-1]     Molecular mass of dry air
    const.MCO2      = 44;         % [g mol-1]     Molecular mass of carbon dioxide
    const.sigmaSB   = 5.67E-8;    % [W m-2 K-4]   Stefan Boltzman constant  
    const.deg2rad   = pi/180;     % [rad]         Conversion from deg to rad
    const.C2K       = 273.15;     % [K]           Melting point of water
    const.deg2rad   = pi/180;     % [rad/deg]     From degrees to radians
end


// src/IO/define_temp_response_biochem.m

function TDP = define_temp_response_biochem

TDP.delHaV	=65330;
TDP.delSV	=485;
TDP.delHdV	=149250;
	
TDP.delHaJ	=43540;
TDP.delSJ	=495;
TDP.delHdJ	=152040;
	
TDP.delHaP	=53100;
TDP.delSP	=490;
TDP.delHdP	=150650;
	
TDP.delHaR	=46390;
TDP.delSR	=490;
TDP.delHdR	=150650;
	
TDP.delHaKc	=79430;
TDP.delHaKo	=36380;
TDP.delHaT	=37830;
	
TDP.Q10	=2;
TDP.s1	=0.3;
TDP.s2	=313.15;
TDP.s3	=0.2;
TDP.s4	=288.15;
TDP.s5	=1.3;
TDP.s6	=328.15;

// src/IO/fill_output_with_nans.m

function [canopy, fluxes, rad, resistance, iter] = fill_output_with_nans(canopy, spectral)

    %% canopy
    
    canopy_names = {"LAIsunlit", "LAIshaded", ...
        "Pntot", "Pnsun", "Pnsha", ...
        "Pntot_Cab", "Pnsun_Cab", "Pnsha_Cab", ...
        "Pntot_Car", "Pnsun_Car", "Pnsha_Car", ...
        "Rntot_PAR", "Rnsun_PAR", "Rnsha_PAR", ...
        "Rntot_Cab", "Rnsun_Cab", "Rnsha_Cab", ...
        "Rntot_Car", "Rnsun_Car", "Rnsha_Car", ...
        "A", "Ja", "ENPQ", "PNPQ", "fqe", "LST", "emis"};
    
    for i=1:length(canopy_names)
        canopy.(canopy_names{i}) = nan;
    end
    
    %% fluxes
    
    fluxes_names = {'Rnctot','lEctot','Hctot','Actot','Tcave', ...
        'Rnstot','lEstot','Hstot','Gtot','Tsave',...
        'Rntot','lEtot','Htot'};
    
    fluxes = struct();
    for i=1:length(fluxes_names)
        fluxes.(fluxes_names{i}) = nan;
    end
    
    %% rad scalars
    
    rad_names = {"PAR", "EPAR",...
        "Eouto","Eoutt","Eoutte", "Lo","Lot","Lote",...
        "F685","wl685","F740","wl740","F684","F761", "LoutF","EoutF","EoutFrc"};
    
    rad = struct();
    for i=1:length(rad_names)
        rad.(rad_names{i}) = nan;
    end
    
    %% rad spectral
    
    rad_names_spectral = {"reflapp", "refl", "rsd", "rdd", "rso", "rdo", ...
        "Eout_", "Lotot_", "Lototf_", "Esun_", "Esky_"};
    
    for i=1:length(rad_names_spectral)
        rad.(rad_names_spectral{i}) = nan(size(spectral.wlS));
    end
    
    %% sif spectral
    
    sif_names = {"LoF_", "sigmaF", "EoutFrc_", "Femleaves_", "EoutF_"};
    
    for i=1:length(sif_names)
        rad.(sif_names{i}) = nan(size(spectral.wlF));
    end
    
    %% resistances
    
    resistances_names = {'raa','raws','rss','ustar'};
    
    resistance = struct();
    for i=1:length(resistances_names)
        resistance.(resistances_names{i}) = nan;
    end

    %% iter
    % in case NaN is in the first row

    iter.counter = nan;

end

// src/IO/initialize_output_structures.m

function [rad,thermal,fluxes] = initialize_output_structures(spectral)

[    fluxes.Rntot,     fluxes.lEtot,      fluxes.Htot,   fluxes.Atot            ,...
    fluxes.Rnctot,    fluxes.lEctot,     fluxes.Hctot,  fluxes.Actot           ,...
    fluxes.Rnstot,    fluxes.lEstot,     fluxes.Hstot,  fluxes.Gtot,   fluxes.Resp    ,...
    thermal.Tcave,    thermal.Tsave   ,...
    thermal.raa,      thermal.rawc,       thermal.raws,   thermal.ustar           ,...
    rad.Lout,      rad.Loutt ,  rad.Eoutte, rad.PAR  ]      =   deal(NaN);
thermal.Ts                             =        NaN(2,1);
%Fc                                 =   deal(NaN(nl,1));

[rad.LoF_      ,...
    rad.Fhem_]           = deal(NaN(size(spectral.wlF,1),1));

[rad.Eouto, rad.Eout ] = deal(NaN);
[rad.Lout_,rad.Lo_]                   =   deal(NaN(size(spectral.wlS,1)),1);
thermal.Ta = NaN;

// src/IO/input_mSCOPE.m

function mly=input_mSCOPE(parameter_file)
    mlayer        = xlsread(parameter_file);
    mly.nly         = mlayer(1,1);
    mly.pLAI        = mlayer(2,:);
    mly.pCab        = mlayer(3,:);
    mly.pCca        = mlayer(4,:);
    mly.pCdm        = mlayer(5,:);
    mly.pCw         = mlayer(6,:);
    mly.pCs         = mlayer(7,:);
    mly.pN          = mlayer(8,:);
    totLAI          = sum(mly.pLAI);
    mly.totLAI      = totLAI;
%     else
%         vi = ones(length(V),1);
%         [~,leafbio,canopy,~,~,~]  = select_input(V,vi,canopy,options,constants);
%         mly.pLAI = canopy.LAI;
%         mly.totLAI = canopy.LAI; 
%         mly.mCab = leafbio.Cab;
%         mly.mCca = leafbio.Cca;
%         mly.mCdm = leafbio.Cdm;
%         mly.mCw = leafbio.Cw;
%         mly.mCs = leafbio.Cs;
%         mly.mN = leafbio.N;
%     end

%% just in case averages
%     mly.mCab        = (mly.pCab*mly.pLAI')./totLAI;
%     mly.mCca        = (mly.pCca*mly.pLAI')./totLAI;
%     mly.mCdm        = (mly.pCdm*mly.pLAI')./totLAI;
%     mly.mCw         = (mly.pCw*mly.pLAI')./totLAI;
%     mly.mCs         = (mly.pCs*mly.pLAI')./totLAI;
%     mly.mN          = (mly.pN*mly.pLAI')./totLAI;
end

// src/IO/load_atmo.m

function atmo = load_atmo(atmfile, SCOPEspec)
    % default error is not clear enough
    assert(exist(atmfile, 'file') == 2, 'Atmospheric file `%s` does not exist', atmfile)
    [~, ~, ext] = fileparts(atmfile);
    if strcmp(ext, '.atm')
        atmo.M  = aggreg(atmfile, SCOPEspec);
    else
        raddata = load(atmfile);
        atmo.Esun_ = raddata(:,1);
        atmo.Esky_ = raddata(:,2);
    end
end

// src/IO/load_mSCOPE_ts.m

function mly_ts = load_mSCOPE_ts(mSCOPE_ts_path, nly, t_column, t_)

%     mly_df = readtable("D:\PyCharm_projects\SCOPE2.0\input\dataset LatinHypercube\input_mSCOPE_timeseries.csv");
    mly_df = readtable(mSCOPE_ts_path);

    t_mly = mly_df.(t_column);
    if all(t_mly <= 367)  % doy is provided
        year_n = unique(year(t_));
        assert(length(year_n) == 1, 'Multiple seasons in mSCOPE are supported only if t_column is TIMESTAMP, not DOY')
        t_mly = datestr(datenum(year_n, 0, t_mly), 'yyyymmddHHMMSS.FFF');
    end
    t_mly = timestamp2datetime(t_mly);

    assert(min(t_mly) <= min(t_) && max(t_mly) >= max(t_), ['Interpolation of mSCOPE values is not possible, '...
        'time range from mSCOPE file does not fully cover range of TS data'])
    
    mly_df.(t_column) = [];
    variables = mly_df.Properties.VariableNames;
    param_names = variables(1:nly:length(variables));

    mly_in_t_ = interp1(t_mly, table2array(mly_df), t_);

    % length(t_) == size(mly_in_t_, 1)
    mly_split = mat2cell(mly_in_t_, length(t_), repmat(nly, 1, length(param_names)));

    mly_ts = cell2struct(mly_split, param_names, 2);
end

// src/IO/load_timeseries.m

function [V, xyt, mly_ts, atmo_paths]  = load_timeseries(V, F, xyt, path_input)
    
    %% filenames
    Dataset_dir = ['dataset ' F(5).FileName];
    meteo_ec_csv = F(6).FileName;
    vegetation_retrieved_csv  = F(7).FileName;
    mSCOPE_csv = F(10).FileName;

    t_column = F(strcmp({F.FileID}, 't')).FileName;
    year_column = F(strcmp({F.FileID}, 'year')).FileName;
    atmo_column = F(strcmp({F.FileID}, 'atmos_names')).FileName;
    
    %% read berkeley format dataset
    df = readtable(fullfile(path_input, Dataset_dir, meteo_ec_csv), ...
        'TreatAsEmpty', {'.','NA','N/A'});
%     df = standardizeMissing(df, -9999); > 2013a
    t_ = df.(t_column);

    if iscell(t_)
        t_ = cellfun(@str2num, t_);
    end
    t_(t_ == -9999) = nan;
    
    if all(t_ <= 367)  % doy is provided
        %assert(~isempty(year_column), 'Please, provide year in your .csv')
        if(isempty(year_column)) 
            warning('t is DOY, converting to date with year = 2020, as `year` in .csv was empty')
            year_n = 2020; 
        else
            year_n = df.(year_column);
        end
        % then we calculate ts for you
        t_ = datestr(datenum(year_n, 0, t_), 'yyyymmddHHMMSS.FFF');
    end
    
    t_ = timestamp2datetime(t_);

    xyt.startDOY = timestamp2datetime(xyt.startDOY);
    xyt.endDOY = timestamp2datetime(xyt.endDOY);
    year_n = year(t_);
    
    %% filtering

    assert(~any(isnat(t_)), 'some t_ in timeseries contain NaN, please, remove or provide. You must know the date')
    
    % this filtering also removes NaNs in t_
    % may be not desirable if df_ts is going to be merged directly => assert
    time_i = (t_ >= xyt.startDOY) & (t_ <= xyt.endDOY);
    df_sub = df(time_i, :);

    %% time 
    t_ = t_(time_i);
    xyt.t = t_;
    xyt.year = year_n(time_i);  % for legacy and doy to date convertion

    %% optional interpolation_csv file
    interpolatable_cols = {};
    if ~isempty(vegetation_retrieved_csv)
        df_int = readtable(fullfile(path_input, Dataset_dir, vegetation_retrieved_csv), ...
            'TreatAsEmpty', {'.','NA','N/A'});
        t_int = df_int.(t_column);
        t_int = timestamp2datetime(t_int, year_n);
        if (min(t_int) > max(t_)) || (max(t_int) < min(t_))
            error('Correct the timestamp, please. The timestamp of vegetation_retrieved_csv is outside of the timestamp of meteo_ec_csv, all simulations would be 0 or NaNs')
        elseif (min(t_int) > min(t_)) && (max(t_int) < max(t_))
            warning('the timestamp of vegetation_retrieved_csv starts later and ends earlier than the timestamp of meteo_ec_csv, head and tail simulations will be 0 or NaNs')
        elseif min(t_int) > min(t_)
            warning('the timestamp of vegetation_retrieved_csv starts later than the timestamp of meteo_ec_csv, head simulations will be 0 or NaNs')
        elseif max(t_int) < max(t_)
            warning('the timestamp of vegetation_retrieved_csv ends earlier than the timestamp of meteo_ec_csv, tail simulations will be 0 or NaNs')
        end
        interpolatable_cols = df_int.Properties.VariableNames;
    end
    
    %% optional mSCOPE file
    mly_ts = struct();
    if ~isempty(mSCOPE_csv)
        mSCOPE_ts_path = fullfile(path_input, Dataset_dir, mSCOPE_csv);
        nly = str2double(F(11).FileName);
        mly_ts = load_mSCOPE_ts(mSCOPE_ts_path, nly, t_column, t_);
        mly_ts.nly = nly;
    end
    
    %% make correspondence: F.FileID : index in V struct
    i_empty = cellfun(@isempty, {F.FileName});
    f_ids = {F(~i_empty).FileID};
    f_names = {F(~i_empty).FileName};
    v_names = {V.Name};
    [~, iF, iV] = intersect(f_ids, v_names, 'stable');
    
    %% read fields that were provided (f_ids)
    for i = 1:length(iF)  % == length(iV)
        fi_i = iF(i);
        vi_i = iV(i);
        col_name = char(f_names(fi_i));
        % TODO replace missing by default?
        if any(strcmp(interpolatable_cols, col_name))
            V(vi_i).Val = interp1(t_int, df_int.(col_name), t_);
        else
            tmp = df_sub.(col_name);
            tmp(tmp == -9999) = nan;
            if all(isnan(tmp))
                warning('%s has NaNs along all timestamp. Calculations may fail', col_name)
            end
            V(vi_i).Val = tmp;
        end
    end

    %% special cases
    %% Irradiance files
    atmo_paths = [];
    if ~isempty(atmo_column)
        atmo_paths = fullfile(path_input, Dataset_dir, df_sub.(atmo_column));
    end
    %% tts calculation
    if ~any(strcmp(f_ids, 'tts'))  % tts wasn't read
        vi_tts = strcmp(v_names, 'tts');
        if isdatetime(t_)
            get_doy = @(x) juliandate(x) - juliandate(datetime(year(x), 1, 0));
            t_ = get_doy(t_);
        end
        DOY_  = floor(t_);
        time_ = 24*(t_-DOY_);
        if all(time_ == 0)
            time_ = 12;
			      warning('taking midday time to calculate tts')
        end
        ttsR  = calczenithangle(DOY_,time_ - xyt.timezn ,0,0,xyt.LON,xyt.LAT);     %sun zenith angle in rad
        V(vi_tts).Val = min(85, ttsR / pi * 180);     
    end

    %% ea calculation
    if ~any(strcmp(f_ids, 'ea')) && any(strcmp(f_ids, 'Ta'))  % ea wasn't read but Ta was
        ta = V(strcmp(v_names, 'Ta')).Val;
        es = satvap(ta);
        vi_ea = strcmp(v_names, 'ea');
        rh_column = F(strcmp({F.FileID}, 'RH')).FileName;
        vpd_column = F(strcmp({F.FileID}, 'VPD')).FileName;
        if ~isempty(rh_column)
            rh = df_sub.(rh_column);
            rh(rh == -9999) = nan;
            if any(rh > 10)
                rh = rh / 100;    % rh from [0 100] to [0 1]
                warning('converted relative hudimity from [0 100] to [0 1]')
            end
            ea = es .* rh;
            warning('calculated ea from Ta and RH')
            V(vi_ea).Val = ea;
        elseif ~isempty(vpd_column)
            vpd = df_sub.(vpd_column);
            vpd(vpd == -9999) = nan;
            ea = es - vpd;
            warning('calculated ea from Ta and VPD')
            if any(ea < 0)
                warning('some ea < 0, is your VPD in hPa?')
            end
            V(vi_ea).Val = ea;
        end 
    end

    %% units convertion
    %% p
    if any(strcmp(f_ids, 'p'))
        vi_p = strcmp(v_names, 'p');
        p = V(vi_p).Val;
        if any(p < 500)
            p = p * 10;
            warning('converted air pressure from kPa to hPa')
        end
        V(vi_p).Val = p;
    end

    %% smc
    if any(strcmp(f_ids, 'SMC'))
        vi_smc = strcmp(v_names, 'SMC');
        smc = V(vi_smc).Val;
        if any(smc > 1)
            smc = smc / 100;  % SMC from [0 100] to [0 1]
            warning('converted soil moisture content from from [0 100] to [0 1]')
        end     
        V(vi_smc).Val = smc;
    end
end

// src/IO/output_data_binary.m

function n_col = output_data_binary(f, k, xyt, rad,  canopy, V, vi, vmax, options, fluxes, meteo, iter,resistance)
%% OUTPUT DATA
% author C. Van der Tol
% date:      30 Nov 2019
% update:    22 Jan 2021 (additional output variables)

%%
if isdatetime(xyt.t)
    get_doy = @(x) juliandate(x) - juliandate(datetime(year(x), 1, 0));
    V(46).Val = get_doy(timestamp2datetime(xyt.startDOY));
    V(47).Val = get_doy(timestamp2datetime(xyt.endDOY));
    xyt.t = get_doy(xyt.t);
end

%% Vegetation products
apar_out= [k xyt.year(k) xyt.t(k)  rad.PAR rad.EPAR canopy.LAIsunlit  canopy.LAIshaded...
    canopy.Pntot     canopy.Pnsun     canopy.Pnsha ...
    canopy.Pntot_Cab canopy.Pnsun_Cab canopy.Pnsha_Cab ...
    canopy.Pntot_Car canopy.Pnsun_Car canopy.Pnsha_Car ...
    canopy.Rntot_PAR canopy.Rnsun_PAR canopy.Rnsha_PAR ...
    canopy.Rntot_Cab canopy.Rnsun_Cab canopy.Rnsha_Cab ...
    canopy.Rntot_Car canopy.Rnsun_Car canopy.Rnsha_Car];
n_col.apar = length(apar_out);
fwrite(f.apar_file,apar_out,'double');

veg_out = [k xyt.year(k) xyt.t(k) canopy.A canopy.Ja canopy.ENPQ  canopy.PNPQ canopy.fqe canopy.LST canopy.emis canopy.GPP];
n_col.veg = length(veg_out);
fwrite(f.veg_file,veg_out,'double');

%% Fluxes product
flu_out = [k iter.counter xyt.year(k) xyt.t(k) cell2mat(struct2cell(fluxes))'];
n_col.flu = length(flu_out);
fwrite(f.flu_file,flu_out,'double');

%% Radiation
rad_out = [k xyt.year(k) xyt.t(k) meteo.Rin, meteo.Rli, rad.Eouto, rad.Eoutt + rad.Eoutte, ...
    rad.Lo, rad.Lot, rad.Lote];
n_col.rad = length(rad_out);
fwrite(f.rad_file,rad_out,'double');

%% Fluorescence scalar outputs
if options.calc_fluor
    fluor_out = [rad.F685  rad.wl685 rad.F740 rad.wl740 rad.F684 rad.F761 ...
        rad.LoutF rad.EoutF rad.EoutFrc];
    n_col.fluor = length(fluor_out);
    fwrite(f.fluor_file,fluor_out,'double');

    %% Fluorescence spectral outputs
    % fluorescence radiance (L) in observation direction [mW m-2 nm-1 sr-1]
    n_col.fluor_spectrum = length(rad.LoF_);
    fwrite(f.fluor_spectrum_file, rad.LoF_, 'double');

    n_col.sigmaF = length(rad.sigmaF);
    fwrite(f.sigmaF_file, rad.sigmaF, 'double');

    n_col.fRC = length(rad.EoutFrc_);
    fwrite(f.fRC_file, rad.EoutFrc_, 'double');

    n_col.fRCL = length(rad.Femleaves_);
    fwrite(f.fRCL_file, rad.Femleaves_, 'double');

    n_col.fhemis = length(rad.EoutF_);
    fwrite(f.fhemis_file,rad.EoutF_, 'double');

    n_col.Lo2 = length(rad.Lototf_);
    fwrite(f.Lo2_file, rad.Lototf_,'double');

    n_col.rapp = length(rad.reflapp);
    fwrite(f.rapp_file, rad.reflapp,'double');


end

%% reflectance
if isfield(f, 'r_file')
    n_col.r = length(rad.refl);
    fwrite(f.r_file,rad.refl,'double');

    n_col.rsd = length(rad.rsd);
    fwrite(f.rsd_file,rad.rsd,'double');

    n_col.rdd = length(rad.rdd);
    fwrite(f.rdd_file,rad.rdd,'double');

    n_col.rso = length(rad.rso);
    fwrite(f.rso_file,rad.rso,'double');

    n_col.rdo = length(rad.rdo);
    fwrite(f.rdo_file,rad.rdo,'double');

    %% Radiance
    n_col.Eout = length(rad.Eout_);
    fwrite(f.Eout_file,rad.Eout_,'double');

    n_col.Lo = length(rad.Lotot_);
    fwrite(f.Lo_file, rad.Lotot_, 'double');

    n_col.Esun = length(rad.Esun_);
    fwrite(f.Esun_file, rad.Esun_, 'double');

    n_col.Esky = length(rad.Esky_);
    fwrite(f.Esky_file, rad.Esky_, 'double');
end

%% Resistances
resist_out = [resistance.raa, resistance.raws, resistance.rss, resistance.ustar];
n_col.resist = length(resist_out);
fwrite(f.resist_file, resist_out, 'double');

%% pars
k2 = find(vmax>1);  % really needed for the first one, later vi > 1
V_short = nan(1,length(k2)+1);
V_short(1) = length(k2);
for i = 1:length(k2)
    V_short(i+1) = V(k2(i)).Val(vi(k2(i)));
end
n_col.pars = length(V_short);
fwrite(f.pars_file, V_short,'double');

end

// src/IO/output_verification_csv.m

function output_verification_csv(Output_dir, verification_dir)
% Date: 07 August 2012
% Author: Christiaan van der Tol (c.vandertol@utwente.nl)
% output_verification.m (script) checks if the output of the latest run
% with SCOPE_v1.51 matches with a 'standard' output located in a directory
% called 'verificationdata'. If it does not, warnings will appear in the
% Matlab command window.
% The following is tested:
%   - does the number of output files match?
%   - does the size of the files match (number of bytes)?
%   - are all files that are in the verification dataset present with the
%   same file names?
%   - is the content of the files exactly the same?
% If the output is different, for example because different parameter
% values have been used in the simulations, then the variables that are
% different will be plotted: the verification data in blue, and the latest
% run in red. In this way the differences can be visually inspected.

% clc, close all
%
% directories         =   dir(['..' filesep 'output' filesep '*']);
% [time_value_s,I]    =   sort([directories(3:end).datenum]);
% Directory           =   directories(2+I(end)).name;
%
% Directory = Output_dir

%% load verification data
path0_ = ['output' filesep verification_dir filesep];
path1_ = [Output_dir filesep];

info0   = dir([path0_ filesep '*.csv']);         %'standard' validation data (to compare with)
info1   = dir([path1_ filesep '*.csv']);           %the most recent output

[differentsize,differentcontent,differentnumberoffiles]  = deal(0);

if ~(length(info0)==length(info1))
    fprintf(['\nWarning: in the output file, ' num2str(length(info1)) ' files were stored, \r'])
    fprintf(['whereas there should be ' num2str(length(info0)) ' files in this directory \r '])
    fprintf('check the simulation options that are specified the options tab of the input spreadsheet \r')
    differentnumberoffiles = 1;
end

L = length(info0);
for i = 1:L
    s0 = info0(i).bytes;
    n0 = info0(i).name;
    for j = 1:length(info1)
        k = strcmp(info1(j).name, n0);
        if k, break, end
    end
    if k
        s1 = info1(j).bytes;
        if ~(s0==s1)
            fprintf(['\n Warning: the file size of ' info0(i).name ' is different from the verification output \r'])
            fprintf(['(' num2str(s1) ' instead of ' num2str(s0) ' bytes) \r'])
            differentsize = 1;
        end
        if (~strcmp(info0(i).name,'pars_and_input.csv') && ~strcmp(info0(i).name,'pars_and_input_short.csv'))
            D0 = dlmread([path0_ info0(i).name],',',2,0);
            D1 = dlmread([path1_ info1(j).name],',',2,0);
        elseif strcmp(info0(i).name,'pars_and_input_short.csv')
            continue
        else
            D0 = dlmread([path0_ info0(i).name],',',1,0);
            D1 = dlmread([path1_ info1(j).name],',',1,0);
        end
        if size(D0) == size(D1)
            if (sum(sum(D0-D1).^2))>1E-9
                fprintf(['\nWarning: data in the output file ' info0(i).name ' are different from the verification output \r '])
                h0 = textread([path0_ info0(i).name],'%s','bufsize', 1E9); %#ok<DTXTRD>
                spn = ceil(sqrt(size(D0,2)));
                figure(i)
                if spn>7
                    nr = min(size(D1, 1), size(D0, 1));
                    for z = 1:nr
                        plot(D0(z,:)','k'), hold on, plot(D1(z,:)','r')
                    end
                    title(info0(i).name, 'interpreter', 'none')
                else
                    h0 = strsplit(h0{1}, ',');
                    for m = 1:size(D0,2)
                        subplot(spn,spn,m)
                        plot(D0(:,m),'k'), hold on, plot(D1(:,m),'r')
                        title([info0(i).name h0(m)], 'interpreter', 'none')
                    end
                end
                differentcontent = 1;
            end
           % differentcontent = 1;
        end
        %end
    else
        fprintf(['\nWarning: the file ' info0(i).name ' was not found in the output\r'])
    end
end

if differentsize
    fprintf('\nWarning The size of some of the output files is different from the verification data \r')
    if differentcontent
        fprintf('Check if the startdate and enddate in the spreadsheet\r')
        fprintf('and the verification data in  are specified in "Dataset_Dir" in the Filenames tab of the input data spreadsheet \r')
    else
        fprintf('but the numerical values are the same. Possible cause: Different Matlab version\r')
    end
end
if ~(differentsize || differentcontent || differentnumberoffiles)
    fprintf('The output is the same as in the verification data set \r')
end
return

// src/IO/savebrdfoutput.m

function savebrdfoutput(options,directional,angles,spectral,Output_dir)
%% directional (This is never saved as binary, but always as ASCII).

if ~exist([Output_dir '\directional'], 'dir')
    mkdir([Output_dir '\directional'])
end
    
Output_angle    =   [directional.tto';  directional.psi']; %#ok<*NASGU>
Output_refl     =   [spectral.wlS    directional.refl_];
Output_rso      =   [spectral.wlS     directional.rso_];
if options.calc_planck
    Output_rad      =   [spectral.wlT'   directional.Lot_];
end
if options.calc_fluor
    Output_fluor = [spectral.wlF'     directional.LoF_];
end

save([Output_dir,'Directional/',sprintf('refl (SunAngle %2.2f degrees).dat',angles.tts)],'Output_refl' ,'-ASCII','-TABS')
save([Output_dir,'Directional/',sprintf('rso (SunAngle %2.2f degrees).dat',angles.tts)],'Output_rso' ,'-ASCII','-TABS')
save([Output_dir,'Directional/',sprintf('Angles (SunAngle %2.2f degrees).dat',angles.tts)],'Output_angle','-ASCII','-TABS')
if options.calc_planck
    save([Output_dir,'Directional/',sprintf('Thermal radiances (SunAngle %2.2f degrees).dat',angles.tts)],'Output_rad','-ASCII','-TABS')
end
if options.calc_fluor
    save([Output_dir,'Directional/',sprintf('Fluorescence (SunAngle %2.2f degrees).dat',angles.tts)],'Output_fluor','-ASCII','-TABS')
end

fiddirtir       =   fopen([Output_dir,'Directional/','read me.txt'],'w');
fprintf(fiddirtir,'The Directional data is written in three files: \r\n');
fprintf(fiddirtir,'\r\n- Angles: contains the observatioin angles. \r\n');
fprintf(fiddirtir,'   * The 1st row gives the observation zenith  angles\r\n');
fprintf(fiddirtir,'   * The 2nd row gives the observation azimuth angles\r\n');
fprintf(fiddirtir,'\r\n- Thermal radiances: contains the radiance at spectral.wlT, emitted plus reflected incoming \r\n');
fprintf(fiddirtir,'   * The 1st column gives the wl values \r\n');
fprintf(fiddirtir,'   * The 2nd column gives the radiances corresponding to the directions given by first column in the Angles file\r\n');
fprintf(fiddirtir,'\r\n- refl and rso: contains the bidirectional distribution functions values, reflectance (based on rso and rdo) and rso. \r\n');
fprintf(fiddirtir,'   * The 1st column gives the wl values corresponding to the BRDF values\r\n');
fprintf(fiddirtir,'   * The 2nd column gives the reflectance values corresponding to the directions given by first column in the Angles file\r\n');
fclose(fiddirtir);

// src/IO/select_input.m

function [soil,leafbio,canopy,meteo,angles,xyt] = select_input(V,vi,canopy,options,constants,xyt,soil,leafbio)

soil.spectrum           = V(60).Val(vi(60));
soil.rss                = V(64).Val(vi(64));
soil.rs_thermal         = V(65).Val(vi(65));
soil.cs                 = V(66).Val(vi(66));
soil.rhos               = V(67).Val(vi(67));
soil.CSSOIL             = V(43).Val(vi(43));
soil.lambdas            = V(68).Val(vi(68));
soil.rbs                = V(44).Val(vi(44));
soil.SMC                = V(54).Val(vi(54));
soil.BSMBrightness      = V(61).Val(vi(61));
soil.BSMlat             = V(62).Val(vi(62));
soil.BSMlon             = V(63).Val(vi(63));

leafbio.Cab             = V(1).Val(vi(1));
leafbio.Cca             = V(2).Val(vi(2));
if options.Cca_function_of_Cab
   leafbio.Cca = 0.25*V(1).Val(vi(1));
end
leafbio.Cdm             = V(3).Val(vi(3));
leafbio.Cw              = V(4).Val(vi(4));
leafbio.Cs              = V(5).Val(vi(5));
leafbio.N               = V(6).Val(vi(6));
leafbio.Cant            = V(7).Val(vi(7));
leafbio.Cp              = V(8).Val(vi(8));
leafbio.Cbc             = V(9).Val(vi(9));
leafbio.rho_thermal     = V(10).Val(vi(10));
leafbio.tau_thermal     = V(11).Val(vi(11));

leafbio.Vcmax25         = V(12).Val(vi(12));
leafbio.BallBerrySlope  = V(13).Val(vi(13));
leafbio.BallBerry0      = V(14).Val(vi(14)); 
leafbio.Type            = V(15).Val(vi(15));
canopy.kV               = V(16).Val(vi(16));
leafbio.RdPerVcmax25    = V(17).Val(vi(17));
fqe                     = V(18).Val(vi(18));
leafbio.Kn0             = V(19).Val(vi(19));
leafbio.Knalpha         = V(20).Val(vi(20));
leafbio.Knbeta          = V(21).Val(vi(21));

leafbio.Tyear           = V(55).Val(vi(55));
leafbio.beta            = V(56).Val(vi(56));
leafbio.kNPQs           = V(57).Val(vi(57));
leafbio.qLs             = V(58).Val(vi(58));
leafbio.stressfactor    = V(59).Val(vi(59));

canopy.LAI              = max(1E-9,V(22).Val(vi(22)));
canopy.hc               = V(23).Val(vi(23));
canopy.zo               = V(24).Val(vi(24));
canopy.d                = V(25).Val(vi(25));
canopy.LIDFa            = V(26).Val(vi(26));
canopy.LIDFb            = V(27).Val(vi(26)); % this is correct (26 instead of 27)
canopy.leafwidth        = V(28).Val(vi(28));
canopy.rb               = V(38).Val(vi(38));
canopy.Cd               = V(39).Val(vi(39));
canopy.CR               = V(40).Val(vi(40));
canopy.CD1              = V(41).Val(vi(41));
canopy.Psicor           = V(42).Val(vi(42));
canopy.rwc              = V(45).Val(vi(45));

%canopy.Cv = V(65).Val(vi(65));
%canopy.crowndiameter = V(66).Val(vi(66));

meteo.z             = V(29).Val(vi(29));
meteo.Rin           = V(30).Val(vi(30));
meteo.Ta            = V(31).Val(vi(31));
meteo.Rli           = V(32).Val(vi(32));
meteo.p             = V(33).Val(vi(33));
meteo.ea            = V(34).Val(vi(34));
meteo.u             = V(35).Val(vi(35));
meteo.Ca            = V(36).Val(vi(36));
meteo.Oa            = V(37).Val(vi(37));

xyt.startDOY        = V(46).Val(vi(46));
xyt.endDOY          = V(47).Val(vi(47));
xyt.LAT             = V(48).Val(vi(48));
xyt.LON             = V(49).Val(vi(49));
xyt.timezn          = V(50).Val(vi(50));

angles.tts          = V(51).Val(vi(51));
angles.tto          = V(52).Val(vi(52));
angles.psi          = V(53).Val(vi(53));

if mean(soil.SMC) > 1
    soil.SMC        = soil.SMC / 100;  % SMC from [0 100] to [0 1]
end   

%% derived input
if options.soil_heat_method ==1
    soil.GAM =  Soil_Inertia1(soil.SMC);
else
    soil.GAM  = Soil_Inertia0(soil.cs,soil.rhos,soil.lambdas);
end

if options.calc_rss_rbs
    [soil.rss,soil.rbs] = calc_rssrbs(soil.SMC,canopy.LAI,soil.rbs);
end

if leafbio.Type
    leafbio.Type = 'C4';
else
    leafbio.Type = 'C3';
end
canopy.hot  = canopy.leafwidth/canopy.hc;
[canopy.zo,canopy.d ]  = zo_and_d(soil,canopy,constants);
leafbio.fqe = fqe;
end

// src/IO/timestamp2datetime.m

function dt = timestamp2datetime(ts, year_n)
    % ts : char or char array like ['20190101'; '20190102']
    %      int or int array [20190101; 20190102]
    % dt : datetime or datetime array

    if iscell(ts)
        ts = cellfun(@str2num, ts);
    end
    ts(ts == -9999) = nan;
    
    if isnumeric(ts) && all(ts <= 367)  % doy is provided
        warning('t is DOY, converting to date with year = %d, as `year` in .csv was empty', year_n(1))
        ts = datestr(datenum(year_n, 0, ts), 'yyyymmddHHMMSS.FFF');
    end
    
    if isnumeric(ts)
        ts = num2str(ts);
        % char or char array like ['20190101'; '20190102']
    end
    
    switch size(ts, 2)
        case 18  % SS.SSS
%             error('milliseconds timestamp is not implemented due to ambiguity of reading: incomplete ts silently become NaNs')
            format = 'yyyyMMddHHmmss.SSS';
        case 14  % SS
            format = 'yyyyMMddHHmmss';
        case 12  % HH, HR
            format = 'yyyyMMddHHmm';
        case 10  % HH (not standard)
            format = 'yyyyMMddHH';
        case 8  % DD
            % in this case HH == 0, MM == 0, SS == 0
            format = 'yyyyMMdd';
        case 6 % MM, WW
            error('weekly and monthly timestamps are not implemented')
        case 4 % YY
            error('yearly timestamp is not implemented')
        otherwise
            error('format of timestamp is not any truncation of `yyyyMMddHHmm`')
    end
    dt = datetime(ts, 'InputFormat', format);
end

// src/RTMs


function rwet = BSM(soilpar,spec,emp)
    
    % Spectral parameters
    
    %wl  = spec.wl;          % wavelengths
    GSV = spec.GSV;         % Global Soil Vectors spectra (nwl * 3)
    kw  = spec.Kw;          % water absorption spectrum
    nw  = spec.nw;          % water refraction index spectrum
    
    % Soil parameters
    
    B   = soilpar.BSMBrightness;        % soil brightness
    lat = soilpar.BSMlat;      % spectral shape latitude (range = 20 - 40 deg)
    lon = soilpar.BSMlon;      % spectral shape longitude (range = 45 - 65 deg)
    SMp = soilpar.SMC*1E2;      % soil moisture volume percentage (5 - 55)
    
    % Empirical parameters
    
    SMC  = emp.SMC;         % soil moisture capacity parameter
    film = emp.film;        % single water film optical thickness
    
    f1 = B * sind(lat);
    f2 = B * cosd(lat) * sind(lon);
    f3 = B * cosd(lat) * cosd(lon);
    
    rdry = f1 * GSV(:,1) + f2 * GSV(:,2) + f3 * GSV(:,3);
    
    % Soil moisture effect
    
    rwet = soilwat(rdry,nw,kw,SMp,SMC,film);
    
return 

function rwet = soilwat(rdry,nw,kw,SMp,SMC,deleff)

    % In this model it is assumed that the water film area is built up  
    % according to a Poisson process. The fractional areas are as follows:
    
    % P(0)   = dry soil area
    % P(1)   = single water film area
    % P(2)   = double water film area
    % ...
    % et cetera
    
    % The fractional areas are given by P(k) = mu^k * exp(-mu) / k! 
    
    % For water films of multiple thickness only the transmission loss due
    % to water absorption is modified, since surface reflectance effects 
    % are not influenced by the thickness of the film
    
    % Input parameters:
    
    % rdry   = dry soil reflectance                             [NW,1]
    % nw     = refraction index of water                        [NW,1]
    % kw     = absorption coefficient of water                  [NW,1]
    % SMp    = soil moisture volume percentage                  [1,NS]
    % SMC    = soil moisture capacity (recommended 0.25)        [1,1]
    % deleff = effective optical thickness of single water film [1,1]
    %          (recommended 0.015)
    
    % Output
    
    % rwet   = wet soil spectra                                 [NW,NS]
    
    % If SMp is given as a row-vector and rdry, nw, kw as column vectors 
    % of the same size (NW, # of wavelengths), then the output is a matrix 
    % of spectra for the different SMp, where each column is a spectrum
    
    % Wout Verhoef
    % Version 1.0
    % September 2012
    
% Peiqi Yang
% Version 1.1
% Jan 2020
% Note by PY (p.yang@utwente.nl)
% because mu<=2, P(k>6) is negligible

%---------------------------------------------------------------------%
k       =   0:6;                    % number of water film, '0' refers to dry soil
nk      =   length(k);              % the number of occurrences
mu      =   (SMp - 5)/ SMC;         % Mu-parameter of Poisson distribution
if mu   <=  0                       % the reason for adding this: if mu<0, fry>1.
    rwet = rdry;                    % we need to check SMC in other parts of SCOPE. soil fluxes routine.
else
    
    % Lekner & Dorf (1988) modified soil background reflectance
    % for soil refraction index = 2.0; uses the tav-function of PROSPECT
    rbac = 1 - (1-rdry) .* (rdry .* tav(90,2.0./nw) / tav(90,2.0) + 1-rdry); % Rbac
    
    % total reflectance at bottom of water film surface
    p    = 1 - tav(90,nw) ./ nw.^2;   % rho21, water to air, diffuse
    
    % reflectance of water film top surface, use 40 degrees incidence angle,
    % like in PROSPECT
    Rw  = 1 - tav(40,nw);             % rho12, air to water, direct
  
    
    % fraction of areas
    % P(0)   = dry soil area            fmul(1)
    % P(1)   = single water film area   fmul(2)
    % P(2)   = double water film area   fmul(3)
    % without loop
    fmul    =   poisspdf(k,mu)';                          % Pobability 
    tw      =   exp(-2*kw * deleff.*k);                   % two-way transmittance,exp(-2*kw*k Delta)
    Rwet_k  =   Rw + (1-Rw) .* (1-p) .*tw.* rbac ./(1 - p .*tw.* rbac);
    rwet   =   rdry * fmul(1) + Rwet_k(:,2:nk)*fmul(2:nk);
end     
return

function Tav = tav(alfa,nr)
n2                                  =   nr.^2;
np                                  =   n2 + 1;
nm                                  =   n2 - 1;

% Stern's formula in Lekner & Dorf (1988) gives reflectance for alfa = 90 degrees

% y1 = (3*n2+2*nr+1)./(3*(nr+1).^2);
% y2 = 2*nr.^3.*(nr.^2+2*nr-1)./(np.^2.*nm);
% y3 = n2.*np.*log(nr)./nm.^2;
% y4 = n2.*nm.^2./np.^3.*log((nr.*(nr+1)./(nr-1)));

% st = y1-y2+y3-y4;

a                                   =   +((nr+1).^2)/2;
k                                   =   -((n2-1).^2)/4;
sin_a                               =   sind(alfa);
%
if alfa~=0    
    B2                              =   sin_a^2 - np/2;
    B1                              =   (alfa~=90) * sqrt( B2.^2 + k );
   
    b                               =   B1 - B2;
    b3                              =   b.^3;
    a3                              =   a.^3;
    
    ts                              =   (k.^2./(6*b3) + k./b - b./2) - ...
                                        (k.^2./(6*a3) + k./a - a./2);
                                    
    tp1                             =   -2*n2.*    (   b  -  a   ) ./ (np.^2);
    tp2                             =   -2*n2.*np.*(  log(b./a)  ) ./ (nm.^2);
    tp3                             =      n2.*    ( 1./b - 1./a ) ./ 2; 
    
%     tp4                             =   16*n2.^2.* (n2.^2+1) .* ( log(2*np.*b - nm.^2) - log(2*np.*a - nm.^2) ) ./ (np.^3.*nm.^2);    
%     tp5                             =   16*n2.^2.* (n2     ) .* ( 1./(2*np.*b - nm.^2) - 1./(2*np.*a - nm.^2)) ./ (np.^3       );

    tp4                             =	16*n2.^2.* (n2.^2+1) .* ( log((2*np.*b - nm.^2)./(2*np.*a - nm.^2))  ) ./(np.^3.*nm.^2);
    tp5                             =   16*n2.^2.* (n2     ) .* ( 1./(2*np.*b - nm.^2) - 1./(2*np.*a - nm.^2)) ./(np.^3);							 
    tp                              =   tp1 + tp2 + tp3 + tp4 + tp5;
    Tav                             =   (ts + tp) ./ (2*sin_a.^2);
else
    Tav                             =   4 *nr/((nr+1)*(nr+1));
end
return

// src/RTMs/RTMf.m

function rad = RTMf(constants,spectral,rad,soil,leafopt,canopy,gap,angles,etau,etah)

% function 'RTMf' calculates the spectrum of fluorescent radiance in the
% observer's direction and also the TOC spectral hemispherical upward Fs flux
%
% Authors:  Wout Verhoef and Christiaan van der Tol (c.vandertol@utwente.nl.nl)
% Date:     12 Dec 2007
% Update:   26 Aug 2008 CvdT        small correction to matrices
%           07 Nov 2008 CvdT        changed layout
% Update:   19 Mar 2009 CvdT        major corrections: lines 95-96,
%                                   101-107, and 119-120.
% Update:    7 Apr 2009 WV & CvdT   major correction: lines 89-90, azimuth
%                                   dependence was not there in previous verions (implicit assumption of
%                                   azimuth(solar-viewing) = 0). This has been corrected
% Update:   May-June 2012 WV & CvdT Add calculation of hemispherical Fs
%                                   fluxes
% Update:   Jan-Feb 2013 WV         Inputs and outputs via structures for
%                                   SCOPE Version 1.40
% Update:   Aug-Oct 2016 PY         Re-write the calculation of emitted SIF
%                                   of each layer. It doesnt use loop at
%                                   all. with the function bsxfun, the
%                                   calculation is much faster
% Update:   Oct 2017-Feb 2018 PY    Re-write the RTM of fluorescence
% Update:   Jan 2020 CvdT           Modified to include 'lite' option,
%                                   mSCOPE representation
% Update:   25 Jun 2020 PY          Po, Ps, Pso. fix the problem we have with the oblique angles above 80 degrees

% Table of contents of the function:
%   0       preparations
%       0.0     globals
%       0.1     initialisations
%       0.2     geometric quantities
%       0.3     solar irradiance factor and ext. in obs dir for all leaf angle/azumith classes
%   1.0     calculation of fluorescence flux in observation direction
%
% Usage: [rad] = RTMf(constants,spectral,rad,soil,leafopt,canopy,gap,angles,etau,etah)
%
% The input and output are structures. These structures are further
% specified in a readme file.
%
% Input:
%   spectral    information about wavelengths and resolutions
%   rad         a large number of radiative fluxes: spectrally distributed
%               and integrated, and canopy radiative transfer coefficients.
%   soil        soil properties
%   leafopt     leaf optical properties
%   canopy      canopy properties (such as LAI and height)
%   gap         probabilities of direct light penetration and viewing
%   angles      viewing and observation angles
%   etau,etah   relative fluorescence emission efficiency for sunlit and
%               shaded leaves, respectively
%
% Output:
%   rad         a large number of radiative fluxes: spectrally distributed
%               and integrated, and canopy radiative transfer coefficients.
%               Here, fluorescence fluxes are added


%% 0.1 initialisations
wlS          = spectral.wlS';       % SCOPE wavelengths, make column vectors
wlF          = (640:4:850)';%spectral.wlF';       % Fluorescence wavelengths
wlE          =   (400:5:750)'; %spectral.wlE';    % Excitation wavelengths
[dummy,iwlfi]    = intersect(wlS,wlE); %#ok<ASGLU>
[dummy,iwlfo]    = intersect(wlS,wlF); %#ok<ASGLU>
nf           = length(iwlfo);
nl           = canopy.nlayers;
LAI          = canopy.LAI;
litab        = canopy.litab;
lazitab      = canopy.lazitab;
lidf         = canopy.lidf;
nlazi        = length(lazitab);         % azumith angle
nlinc        = length(litab);           % inclination
nlori        = nlinc * nlazi;           % total number of leaf orientations
layers       = 1:nl;

Ps           = gap.Ps;
Po           = gap.Po;
Pso          = gap.Pso;

Qs          = Ps(1:end-1);

[MpluEmin   ,...
    MpluEplu   , ...
    MminEmin   , ...
    MminEplu   , ...
    MpluEsun   , ...
    MminEsun]       = deal(zeros(nf,nl));

% for speed-up the calculation only uses wavelength i and wavelength o part of the spectrum
Esunf_             = rad.Esun_(iwlfi);
Eminf_             = rad.Emin_(:,iwlfi)';          % transpose into [nwlfo,nl] matrix
Epluf_             = rad.Eplu_(:,iwlfi)';
iLAI               = LAI/nl;                       % LAI of a layer        [1]

Xdd         = rad.Xdd(:,iwlfo);
rho_dd      = rad.rho_dd(:,iwlfo);
R_dd        = rad.R_dd(:,iwlfo);
tau_dd      = rad.tau_dd(:,iwlfo);
vb          = rad.vb(:,iwlfo);
vf          = rad.vf(:,iwlfo);

%% 0.2 geometric quantities
Mb                = leafopt.Mb;
Mf                = leafopt.Mf;

% geometric factors
deg2rad             = constants.deg2rad;
tto                 = angles.tto;
tts                 = angles.tts;
psi                 = angles.psi;
rs                  = soil.refl(iwlfo,:);           % [nwlfo]     soil reflectance
cos_tto             = cos(tto*deg2rad);             % cos observation zenith angle
sin_tto             = sin(tto*deg2rad);             % sin observation zenith angle

cos_tts             = cos(tts*deg2rad);             % cos solar angle
sin_tts             = sin(tts*deg2rad);             % sin solar angle

cos_ttli            = cos(litab*deg2rad);           % cos leaf inclinaation angles
sin_ttli            = sin(litab*deg2rad);           % sin leaf inclinaation angles
cos_phils           = cos(lazitab*deg2rad);         % cos leaf azimuth angles rel. to sun azi
cos_philo           = cos((lazitab-psi)*deg2rad);   % cos leaf azimuth angles rel. to viewing azi

%% 0.3 geometric factors for all leaf angle/azumith classes
cds                 = cos_ttli*cos_tts*ones(1,36) + sin_ttli*sin_tts*cos_phils;  % [nli,nlazi]
cdo                 = cos_ttli*cos_tto*ones(1,36) + sin_ttli*sin_tto*cos_philo;  % [nli,nlazi]
fs                  = cds/cos_tts;                                               % [nli,nlazi]
absfs               = abs(fs);                                                   % [nli,nlazi]
fo                  = cdo/cos_tto;                                               % [nli,nlazi]
absfo               = abs(fo);                                                   % [nli,nlazi]
fsfo                = fs.*fo;                                                    % [nli,nlazi]
absfsfo             = abs(fsfo);                                                 % [nli,nlazi]
foctl               = fo.*(cos_ttli*ones(1,36));                                 % [nli,nlazi]
fsctl               = fs.*(cos_ttli*ones(1,36));                                 % [nli,nlazi]
ctl2                = cos_ttli.^2*ones(1,36);                                    % [nli,nlazi]

% reshape all the variables
absfs               = reshape(absfs,nlori,1);                                    % [nlori,1]
absfo               = reshape(absfo,nlori,1);                                    % [nlori,1]
fsfo                = reshape(fsfo,nlori,1);                                     % [nlori,1]
absfsfo             = reshape(absfsfo,nlori,1);                                  % [nlori,1]
foctl               = reshape(foctl,nlori,1);                                    % [nlori,1]
fsctl               = reshape(fsctl,nlori,1);                                    % [nlori,1]
ctl2                = reshape(ctl2,nlori,1);                                     % [nlori,1]

%% 1.0 calculation of fluorescence flux in observation direction

% fluorescence matrices and efficiencies
[U,Fmin_,Fplu_] =deal(zeros(nl+1,size(leafopt.Mb,1)));

Mplu = 0.5*(Mb+Mf);    % [nwlfo,nwlfi]
Mmin = 0.5*(Mb-Mf);    % [nwlfo,nwlfi]

% in-products: we convert incoming radiation to a fluorescence spectrum using the matrices.
for j = 1:nl
    ep = constants.A*ephoton(wlF*1E-9,constants);
    MpluEmin(:,j)   = ep.*(Mplu(:,:,j)* e2phot(wlE*1E-9,Eminf_(:,j),constants));          % [nwlfo,nl+1]
    MpluEplu(:,j)   = ep.*(Mplu(:,:,j)* e2phot(wlE*1E-9,Epluf_(:,j),constants));          % [nwlfo,nl+1]
    MminEmin(:,j)   = ep.*(Mmin(:,:,j)* e2phot(wlE*1E-9,Eminf_(:,j),constants));          % [nwlfo,nl+1]
    MminEplu(:,j)   = ep.*(Mmin(:,:,j)* e2phot(wlE*1E-9,Epluf_(:,j),constants));          % [nwlfo,nl+1]
    MpluEsun(:,j)   = ep.*(Mplu(:,:,j)* e2phot(wlE*1E-9,Esunf_,constants));               % integration by inproduct
    MminEsun(:,j)   = ep.*(Mmin(:,:,j)* e2phot(wlE*1E-9,Esunf_,constants));               % integration by inproduct
end

laz= 1/36;

if size(etau,2)<2
    etau = repmat(etau,1,13,36);
    etau = permute(etau,[2 3 1]);
end

etau_lidf = bsxfun(@times,reshape(etau,nlori,nl),repmat(lidf*laz,36,1));     %[nlori,nl]
etah_lidf = bsxfun(@times,repmat(etah,1,nlori)',repmat(lidf*laz,36,1));

wfEs      =  bsxfun(@times,sum(bsxfun(@times,etau_lidf,absfsfo)),MpluEsun) +...
    bsxfun(@times,sum(bsxfun(@times,etau_lidf,fsfo)),MminEsun);
sfEs     =  bsxfun(@times,sum(bsxfun(@times,etau_lidf,absfs)),MpluEsun) -...
    bsxfun(@times,sum(bsxfun(@times,etau_lidf,fsctl)),MminEsun);
sbEs     =  bsxfun(@times,sum(bsxfun(@times,etau_lidf,absfs)),MpluEsun) +...
    bsxfun(@times,sum(bsxfun(@times,etau_lidf,fsctl)),MminEsun);
vfEplu_h  = bsxfun(@times,sum(bsxfun(@times,etah_lidf,absfo)),MpluEplu) -...
    bsxfun(@times,sum(bsxfun(@times,etah_lidf,foctl)),MminEplu);
vfEplu_u  = bsxfun(@times,sum(bsxfun(@times,etau_lidf,absfo)),MpluEplu) -...
    bsxfun(@times,sum(bsxfun(@times,etau_lidf,foctl)),MminEplu);
vbEmin_h  = bsxfun(@times,sum(bsxfun(@times,etah_lidf,absfo)),MpluEmin) +...
    bsxfun(@times,sum(bsxfun(@times,etah_lidf,foctl)),MminEmin);
vbEmin_u  = bsxfun(@times,sum(bsxfun(@times,etau_lidf,absfo)),MpluEmin) +...
    bsxfun(@times,sum(bsxfun(@times,etau_lidf,foctl)),MminEmin);
sigfEmin_h  = bsxfun(@times,sum(etah_lidf),MpluEmin) -...
    bsxfun(@times,sum(bsxfun(@times,etah_lidf,ctl2)),MminEmin);
sigfEmin_u  = bsxfun(@times,sum(etau_lidf),MpluEmin) -...
    bsxfun(@times,sum(bsxfun(@times,etau_lidf,ctl2)),MminEmin);
sigbEmin_h  = bsxfun(@times,sum(etah_lidf),MpluEmin) +...
    bsxfun(@times,sum(bsxfun(@times,etah_lidf,ctl2)),MminEmin);
sigbEmin_u  = bsxfun(@times,sum(etau_lidf),MpluEmin) +...
    bsxfun(@times,sum(bsxfun(@times,etau_lidf,ctl2)),MminEmin);
sigfEplu_h  = bsxfun(@times,sum(etah_lidf),MpluEplu) -...
    bsxfun(@times,sum(bsxfun(@times,etah_lidf,ctl2)),MminEplu);
sigfEplu_u  = bsxfun(@times,sum(etau_lidf),MpluEplu) -...
    bsxfun(@times,sum(bsxfun(@times,etau_lidf,ctl2)),MminEplu);
sigbEplu_h  = bsxfun(@times,sum(etah_lidf),MpluEplu) +...
    bsxfun(@times,sum(bsxfun(@times,etah_lidf,ctl2)),MminEplu);
sigbEplu_u  = bsxfun(@times,sum(etau_lidf),MpluEplu) +...
    bsxfun(@times,sum(bsxfun(@times,etau_lidf,ctl2)),MminEplu);

%   Emitted fluorescence
piLs        =   wfEs+vfEplu_u+vbEmin_u;         % sunlit for each layer
piLd        =   vbEmin_h+vfEplu_h;              % shade leaf for each layer
Fsmin       =   sfEs+sigfEmin_u+sigbEplu_u;     % Eq. 29a for sunlit leaf
Fsplu       =   sbEs+sigbEmin_u+sigfEplu_u;     % Eq. 29b for sunlit leaf
Fdmin       =   sigfEmin_h+sigbEplu_h;          % Eq. 29a for shade leaf
Fdplu       =   sigbEmin_h+sigfEplu_h;          % Eq. 29b for shade leaf
Femmin      =   iLAI*bsxfun(@times,Qs', Fsmin) +iLAI* bsxfun(@times,(1-Qs)',Fdmin);
Femplu      =   iLAI*bsxfun(@times,Qs', Fsplu) +iLAI* bsxfun(@times,(1-Qs)',Fdplu);

for j=nl:-1:1      % from bottom to top
    Y(j,:)  =(rho_dd(j,:).*U(j+1,:)+Femmin(:,j)')./(1-rho_dd(j,:).*R_dd(j+1,:));
    U(j,:) =tau_dd(j,:).*(R_dd(j+1,:).*Y(j,:)+U(j+1,:))+Femplu(:,j)';
end

for j=1:nl          % from top to bottom
    Fmin_(j+1,:)  = Xdd(j,:).*Fmin_(j,:)+Y(j,:);
    Fplu_(j,:)  = R_dd(j,:).*Fmin_(j,:)+U(j,:);
end
piLo1     = iLAI*piLs*Pso(1:nl);
piLo2     = iLAI*piLd*(Po(1:nl)-Pso(1:nl));
piLo3     = iLAI*(vb.*Fmin_(layers,:)  + vf.*Fplu_(layers,:))'*Po(1:nl);
piLo4     = rs .* Fmin_(nl+1,:)' * Po(nl+1);

piLtot      = piLo1 + piLo2 + piLo3 + piLo4;
LoF_        = piLtot/pi;
Fhem_       = Fplu_(1,:)';

method = 'spline';  % M2020a name
%if verLessThan('matlab', '9.8')
%    method = 'splines';
%end

rad.LoF_    = interp1(wlF,LoF_,spectral.wlF',method);
rad.EoutF_   = interp1(wlF,Fhem_,spectral.wlF',method);

rad.LoF_sunlit      = interp1(wlF,piLo1/pi,spectral.wlF',method);
rad.LoF_shaded      = interp1(wlF,piLo2/pi,spectral.wlF',method);
rad.LoF_scattered   = interp1(wlF,piLo3/pi,spectral.wlF',method);
rad.LoF_soil        = interp1(wlF,piLo4/pi,spectral.wlF',method);

rad.EoutF   = 0.001 * Sint(Fhem_,wlF);
rad.LoutF   = 0.001 * Sint(LoF_,wlF);

rad.Femleaves_ = interp1(wlF,sum(Femmin+Femplu,2),spectral.wlF',method);

[rad.F685,iwl685]  = max(rad.LoF_(1:55));
rad.wl685 = spectral.wlF(iwl685);
if iwl685 == 55; [rad.F685,rad.wl685] = deal(NaN); end
[rad.F740,iwl740]  = max(rad.LoF_(70:end));
rad.wl740 = spectral.wlF(iwl740+69);
rad.F684  = rad.LoF_(685-spectral.wlF(1));
rad.F761  = rad.LoF_(762-spectral.wlF(1));
return

// src/RTMs/RTMo.m

function [rad,gap,profiles] = RTMo(spectral,atmo,soil,leafopt,canopy,angles,constants,meteo,options)
% calculates the spectra of hemisperical and directional observed visible 
% and thermal radiation (fluxes E and radiances L), as well as the single 
% and bi-directional gap probabilities
%
% the function does not require any non-standard Matlab functions. No
% changes to the code have to be made to operate the function for a
% particular canopy. All necessary parameters and variables are input or
% global and need to be specified elsewhere.
%
% Authors:      Wout Verhoef            (w.verhoef@utwente.nl) 
%               Christiaan van der Tol  (c.vandertol@utwente.nl)
%               Joris Timmermans        ()
%
% updates:      10 Sep 2007 (CvdT)      - calculation of Rn
%                5 Nov 2007             - included observation direction
%               12 Nov 2007             - included abs. PAR spectrum output
%                                       - improved calculation efficiency
%               13 Nov 2007             - written readme lines
%               11 Feb 2008 (WV&JT)     - changed Volscat
%                           (JT)        - small change in calculation Po,Ps,Pso
%                                        - introduced parameter 'lazitab'
%                                       - changed nomenclature
%                                       - Appendix IV: cosine rule
%               04 Aug 2008 (JT)        - Corrections for Hotspot effect in the probabilities
%               05 Nov 2008 (CvdT)      - Changed layout
%               04 Jan 2011 (JT & CvdT) - Included Pso function (Appendix IV)
%                                       - removed the analytical function (for checking)  
%               02 Oct 2012 (CvdT)      - included incident PAR in output
%
%               Jan/Feb 2013 (WV)       - Major revision towards SCOPE version 1.40:
%                                       - Parameters passed using structures
%                                       - Improved interface with MODTRAN atmospheric data
%                                       - Now also calculates 4-stream
%                                         reflectances rso, rdo, rsd and rdd
%                                         analytically 
%               Apri 2013 (CvT)         - improvements in variable names
%                                           and descriptions
%                  Dec 2019 CvdT        mSCOPE representation, lite option
%
%                                      
% Table of contents of the function
%
%   0.      Preparations
%       0.1     parameters
%       0.2     initialisations
%   1.      Geometric quantities
%       1.1     general geometric quantities
%       1.2     geometric factors associated with extinction and scattering
%       1.3     geometric factors to be used later with rho and tau
%       1.4     solar irradiance factor for all leaf orientations
%       1.5     probabilities Ps, Po, Pso
%   2.      Calculation of upward and downward fluxes
%   3.      Outgoing fluxes, hemispherical and in viewing direction, spectrum
%   4.      Net fluxes, spectral and total, and incoming fluxes
%   A1      functions J1 and J2 (introduced for stable solutions)
%   A2      function volscat
%   A3      function e2phot
%   A4      function Pso
%
% references:
%{1} Verhoef (1998), 'Theory of radiative transfer models applied in
%    optical remote sensing of vegetation canopies'. PhD Thesis Univ. Wageninegn
%{2} Verhoef, W., Jia, L., Xiao, Q. and Su, Z. (2007) Unified optical -
%    thermal four - stream radiative transfer theory for homogeneous 
%    vegetation canopies. IEEE Transactions on geoscience and remote 
%    sensing, 45,6.
%{3} Verhoef (1985), 'Earth Observation Modeling based on Layer Scattering
%    Matrices', Remote sensing of Environment, 17:167-175 
%              
% Usage:
% function [rad,gap,profiles] = RTMo(spectral,atmo,soil,leafopt,canopy,angles,meteo,rad,options)  
%
% The input and output are structures. These structures are further
% specified in a readme file.
%
% Input:
%   spectral    information about wavelengths and resolutions
%   atmo        MODTRAN atmospheric parameters
%   soil        soil properties
%   leafopt     leaf optical properties
%   canopy      canopy properties (such as LAI and height)
%   angles      viewing and observation angles
%   meteo       has the meteorological variables. Is only used to correct
%               the total irradiance if a specific value is provided
%               instead of the usual Modtran output.
%   rad         initialization of the structure of the output 'rad'
%   options     simulation options. Here, the option
%               'calc_vert_profiles' is used, a boolean that tells whether 
%               or not to output data of 60 layers separately.
%
% Output:
%   gap         probabilities of direct light penetration and viewing
%   rad         a large number of radiative fluxes: spectrally distributed 
%               and integrated, and canopy radiative transfer coefficients.
%   profiles    vertical profiles of radiation variables such as absorbed
%               PAR.


%% 0. Preparations
deg2rad = constants.deg2rad;    % degree to rad
wl      = spectral.wlS;         % SCOPE wavelengths as a column-vector
nwl     = length(wl);           %
wlPAR   = spectral.wlPAR;       % PAR wavelength range
minPAR  = min(wlPAR);           % min PAR
maxPAR  = max(wlPAR);           % max PAR
Ipar    = find(wl>=minPAR & wl<=maxPAR); % Indices for PAR wavelenghts within wl
tts     = angles.tts;           % solar zenith angle
tto     = angles.tto;           % observer zenith angle
psi     = angles.psi;           % relative azimuth anglee

nl      = canopy.nlayers;       % number of canopy layers (nl)
litab   = canopy.litab;         % SAIL leaf inclibation angles % leaf inclination angles PY
lazitab = canopy.lazitab;       % leaf azimuth angles relative to the sun
nlazi   = canopy.nlazi;         % number of azimuth angles (36)
LAI     = canopy.LAI;           % leaf area index
lidf    = canopy.lidf;          % leaf inclination distribution function
xl      = canopy.xl;            % all levels except for the top
dx      = 1/nl;

rho = leafopt.refl;
tau = leafopt.tran;
kChlrel = leafopt.kChlrel;
kCarrel = leafopt.kCarrel;
rs      =   soil.refl;          % [nwl,nsoils] soil reflectance spectra
epsc    =   1-rho-tau;          % [nl,nwl]        emissivity of leaves
epss    =   1-rs;              % [nwl]        emissivity of soil
iLAI    =   LAI/nl;             % [1]          LAI of elementary layer

% initializations
Rndif                       = zeros(nl,1);                 % [nl+1]         abs. diffuse rad soil+veg
[Pdif,Pndif,Pndif_Cab,Rndif_Cab,Pndif_Car,Rndif_Car,Rndif_PAR]      = deal(zeros(nl,1));           % [nl]           incident and net PAR veg
%[Es_,Emin_,Eplu_]           = deal(zeros(nl+1,nwl));       % [nl+1,nwl]     direct, up and down diff. rad.
[Rndif_]                    = zeros(nl,nwl);               % [nl,nwl]       abs diff and PAR veg.
[Pndif_,Rndif_PAR_]         = deal(zeros(nl,length(wlPAR)));
[Pndif_Cab_,Rndif_Cab_,Pndif_Car_,Rndif_Car_]         = deal(zeros(nl,length(spectral.IwlP)));
%[Puc,Rnuc,Pnuc,Pnuc_Cab,Rnuc_PAR]    = deal(zeros(nli,nlazi,nl));   %

%% 1. Geometric quantities
% 1.1 general geometric quantities. these variables are scalars
cos_tts     = cos(tts*deg2rad);             %           cos solar       angle
tan_tto     = tan(tto*deg2rad);             %           tan observation angle

cos_tto     = cos(tto*deg2rad);             %           cos observation angle
sin_tts     = sin(tts*deg2rad);             %           sin solar       angle
tan_tts     = tan(tts*deg2rad);             %           tan observation angle

psi         = abs(psi-360*round(psi/360));  %           (to ensure that volscatt is symmetric for psi=90 and psi=270)
dso         = sqrt(tan_tts.^2 + tan_tto.^2 - 2*tan_tts.*tan_tto.*cos(psi*deg2rad));

% 1.2 geometric factors associated with extinction and scattering
[chi_s,chi_o,frho,ftau]=volscat(tts,tto,psi,litab);   % volume scattering

cos_ttlo    = cos(lazitab*deg2rad);         % [1,36]    cos leaf azimuth angles

cos_ttli    = cos(litab*deg2rad);           % [13]      cos leaf angles
sin_ttli    = sin(litab*deg2rad);           % [13]      sinus leaf angles

ksli        = chi_s./cos_tts;               % [13]      p306{1} extinction coefficient in direction of sun        per leaf angle
koli        = chi_o./cos_tto;               % [13]      p307{1} extinction coefficient in direction of observer   per leaf angle

sobli       = frho*pi/(cos_tts*cos_tto);    % [13]      pag 309{1} area scattering coefficient fractions
sofli       = ftau*pi/(cos_tts*cos_tto);    % [13]      pag 309{1}
bfli        = cos_ttli.^2;                  % [13]

%integration over angles (using a vector inproduct) -> scalars
k           = ksli'*lidf;                   %           pag 306{1}    extinction coefficient in direction of sun.
K           = koli'*lidf;                   %           pag 307{1}    extinction coefficient in direction of observer
bf          = bfli'*lidf;                   %
sob         = sobli'*lidf;                  %           weight of specular2directional back    scatter coefficient
sof         = sofli'*lidf;                  %           weight of specular2directional forward scatter coefficient
% 1.3 geometric factors to be used later with rho and tau, f1 f2 of pag 304:
% these variables are scalars
sdb         = 0.5*(k+bf);                   % fs*f1
sdf         = 0.5*(k-bf);                   % fs*f2     weight of specular2diffuse     foward  scatter coefficient
ddb         = 0.5*(1+bf);                   % f1^2+f2^2 weight of diffuse2diffuse      back    scatter coefficient
ddf         = 0.5*(1-bf);                   % 2*f1*f2   weight of diffuse2diffuse      forward scatter coefficient
dob         = 0.5*(K+bf);                   % fo*f1     weight of diffuse2directional  back    scatter coefficient
dof         = 0.5*(K-bf);                   % fo*f2     weight of diffuse2directional  forward scatter coefficient

% 1.4 solar irradiance factor for all leaf orientations
Css          = cos_ttli*cos_tts;             % [nli]     pag 305 modified by Joris
Ss          = sin_ttli*sin_tts;             % [nli]     pag 305 modified by Joris
cos_deltas  = Css*ones(1,nlazi) + Ss*cos_ttlo;  % [nli,nlazi]
fs          = abs(cos_deltas/cos_tts);         % [nli,nlazi] pag 305

%% 2. Calculation of reflectance
% 2.1  reflectance, transmittance factors in a thin layer
% the following are vectors with lenght [nl,nwl]
sigb        = ddb*rho + ddf*tau;            % [nl,nwl]     sigmab, p305{1} diffuse     backscatter scattering coefficient for diffuse  incidence
sigf        = ddf*rho + ddb*tau;            % [nl,nwl]     sigmaf, p305{1} diffuse     forward     scattering coefficient for forward  incidence
sb          = sdb*rho + sdf*tau;            % [nl,nwl]     sb,     p305{1} diffuse     backscatter scattering coefficient for specular incidence
sf          = sdf*rho + sdb*tau;            % [nl,nwl]     sf,     p305{1} diffuse     forward     scattering coefficient for specular incidence
vb          = dob*rho + dof*tau;            % [nl,nwl]     vb,     p305{1} directional backscatter scattering coefficient for diffuse  incidence
vf          = dof*rho + dob*tau;            % [nl,nwl]     vf,     p305{1} directional forward     scattering coefficient for diffuse  incidence
w           = sob*rho + sof*tau;            % [nl,nwl]     w,      p309{1} bidirectional scattering coefficent (directional-directional)
a           = 1-sigf;                       % [nl,nwl]     attenuation

%% 3. Flux calculation

% diffuse fluxes within the vegetation covered part
%iLAI    =   LAI/nl; 
tau_ss = repmat(1-k.*iLAI,nl,1);  %REPLACE when LIDF profile ready.
tau_dd = (1-a.*iLAI);
tau_sd = sf.*iLAI;
rho_sd = sb.*iLAI;
rho_dd = sigb.*iLAI;
[R_sd,R_dd,Xss,Xsd,Xdd] = calc_reflectances(tau_ss,tau_sd,tau_dd,rho_dd,rho_sd,rs,nl,nwl);
rdd     = R_dd(1,:)';
rsd     = R_sd(1,:)';
[Esun_,Esky_] = calcTOCirr(atmo,meteo,rdd,rsd,wl,nwl);

[Emins_,Eplus_] = calc_fluxprofile(Esun_,0*Esky_,rs,Xss,Xsd,Xdd,R_sd,R_dd,nl,nwl);
[Emind_,Eplud_] = calc_fluxprofile(0*Esun_,Esky_,rs,Xss,Xsd,Xdd,R_sd,R_dd,nl,nwl);
Emin_ = Emins_+Emind_;
Eplu_ = Eplus_+Eplud_;
%%
% 1.5 probabilities Ps, Po, Pso
Ps          =   exp(k*xl*LAI) ;                                              % [nl+1]  p154{1} probability of viewing a leaf in solar dir
Po          =   exp(K*xl*LAI) ;
Ps(1:nl)    =   Ps(1:nl) *(1-exp(-k*LAI*dx))/(k*LAI*dx);                                      % Correct Ps/Po for finite dx
Po(1:nl)    =   Po(1:nl) *(1-exp(-K*LAI*dx))/(K*LAI*dx);  % Correct Ps/Po for finite dx
q           =   canopy.hot;
Pso         =   zeros(size(Po));
for j=1:length(xl)
    Pso(j,:)=   quad(@(y)Psofunction(K,k,LAI,q,dso,y),xl(j)-dx,xl(j))/dx; %#ok<FREMO>
end
Pso(Pso>Po)= min([Po(Pso>Po),Ps(Pso>Po)],[],2);    %takes care of rounding error
Pso(Pso>Ps)= min([Po(Pso>Ps),Ps(Pso>Ps)],[],2);    %takes care of rounding error
gap.Pso      = Pso;

%%
% 3.3 outgoing fluxes, hemispherical and in viewing direction, spectrum
% in viewing direction, spectral due to diffuse light

% vegetation contribution
piLocd_     = (sum(vb.*Po(1:nl).*Emind_(1:nl,:)) +...
              sum(vf.*Po(1:nl).*Eplud_(1:nl,:)))'*iLAI;
% soil contribution          
piLosd_     = rs.*(Emind_(end,:)'*Po(end)); 

% in viewing direction, spectral due to direct solar light
% vegetation contribution
piLocu_     = (sum(vb.*Po(1:nl).*Emins_(1:nl,:)) +...
              sum(vf.*Po(1:nl).*Eplus_(1:nl,:))+...
              sum(w.*Pso(1:nl).*Esun_'))'*iLAI;
% soil contribution
piLosu_     = rs.* (Emins_(end,:)'*Po(end) + Esun_*Pso(end)); 

piLod_      = piLocd_ + piLosd_;        % [nwl] piRad in obsdir from Esky
piLou_      = piLocu_ + piLosu_;        % [nwl] piRad in obsdir from Eskun
piLoc_      = piLocu_ + piLocd_;        % [nwl] piRad in obsdir from vegetation
piLos_      = piLosu_ + piLosd_;        % [nwl] piRad in obsdir from soil

piLo_       = piLoc_ + piLos_;          % [nwl] piRad in obsdir

Lo_         = piLo_/pi;                 % [nwl] Rad in obsdir
rso         = piLou_./Esun_;            % [nwl] obsdir reflectance of solar beam
rdo         = piLod_./Esky_;            % [nlw] obsir reflectance of sky irradiance
Refl        = piLo_./(Esky_+Esun_);     % [nwl] 
Refl(Esky_<1E-4) = rso(Esky_<1E-4);     % prevents numerical instability in absorption windows
I            = find(Esky_<2E-4*max(Esky_));
Refl(I)      = rso(I);                  % prevents numerical instability in absorption windows

%% 4. net fluxes, spectral and total, and incoming fluxes
%4.1 incident PAR at the top of canopy, spectral and spectrally integrated
P_          = e2phot(wl(Ipar)*1E-9,(Esun_(Ipar)+Esky_(Ipar)),constants);      %
P           = 0.001 * Sint(P_,wlPAR);                                % mol m-2s-1
EPAR_       = Esun_(Ipar)+Esky_(Ipar);
EPAR        = 0.001 * Sint(EPAR_,wlPAR);
%Psun        = 0.001 * Sint(e2phot(wlPAR*1E-9,Esun_(Ipar),constants),wlPAR);   % Incident solar PAR in PAR units
% Incident and absorbed solar radiation

%4.2 Absorbed radiation
%    absorbed radiation in Wm-2         (Asun)
%    absorbed PAR in mol m-2s-1         (Pnsun)
%    absorbed PAR in Wm-2               (Rnsun_PAR)
%    absorbed PAR by Chl in mol m-2s-1  (Pnsun_Cab)

[Asun,Pnsun,Rnsun_PAR,Pnsun_Cab,Rnsun_Cab,Pnsun_Car,Rnsun_Car]= deal(zeros(nl,1));
for j=1:nl
    Asun(j)        = 0.001 * Sint(Esun_.*epsc(j,:)',wl);                                 % Total absorbed solar radiation
    Pnsun(j)       = 0.001 * Sint(e2phot(wlPAR*1E-9,Esun_(Ipar).*epsc(j,Ipar)',constants),wlPAR);  % Absorbed solar radiation in PAR range in moles m-2 s-1
    Rnsun_Cab(j)   = 0.001 * Sint(Esun_(spectral.IwlP)'.*epsc(j,spectral.IwlP).*kChlrel(j,:),spectral.wlP);
    Rnsun_Car(j)   = 0.001 * Sint(Esun_(spectral.IwlP)'.*epsc(j,spectral.IwlP).*kCarrel(j,:),spectral.wlP);
    Rnsun_PAR(j)   = 0.001 * Sint(Esun_(Ipar)'.*epsc(j,Ipar),wlPAR);
    Pnsun_Cab(j)   = 0.001 * Sint(e2phot(spectral.wlP*1E-9,kChlrel(j,:)'.*Esun_(spectral.IwlP).*epsc(j,spectral.IwlP)',constants),spectral.wlP);
    Pnsun_Car(j)   = 0.001 * Sint(e2phot(spectral.wlP*1E-9,kCarrel(j,:)'.*Esun_(spectral.IwlP).*epsc(j,spectral.IwlP)',constants),spectral.wlP);

end

%4.3 total direct radiation (incident and net) per leaf area (W m-2 leaf)

% total direct radiation (incident and net) per leaf area (W m-2 leaf)
%Pdir        = fs * Psun;                        % [13 x 36]   incident
if options.lite
    fs      = lidf'*mean(fs,2);%
    Rndir       = fs * Asun(j);                        % [13 x 36 x nl]   net
    Pndir       = fs * Pnsun(j);                       % [13 x 36 x nl]   net PAR
    Pndir_Cab   = fs * Pnsun_Cab(j);                   % [13 x 36 x nl]   net PAR Cab
    Rndir_Cab   = fs * Rnsun_Cab(j);                   %   net PAR energy units
    Pndir_Car   = fs * Pnsun_Car(j);                   % [13 x 36 x nl]   net PAR Cab
    Rndir_Car   = fs * Rnsun_Car(j);                   %   net PAR energy units
    Rndir_PAR   = fs * Rnsun_PAR(j);                   % [13 x 36 x nl]
else
    [Rndir,Pndir,Pndir_Cab,Rndir_Cab,Rndir_PAR,Pndir_Car,Rndir_Car]= deal(zeros(13,36,nl));
    for j=1:nl
        Rndir(:,:,j)       = fs * Asun(j);                        % [13 x 36 x nl]   net
        Pndir(:,:,j)       = fs * Pnsun(j);                       % [13 x 36 x nl]   net PAR
        Pndir_Cab(:,:,j)   = fs * Pnsun_Cab(j);                   % [13 x 36 x nl]   net PAR Cab
        Rndir_Cab(:,:,j)   = fs * Rnsun_Cab(j);                   %   net PAR energy units
        Pndir_Car(:,:,j)   = fs * Pnsun_Car(j);                   % [13 x 36 x nl]   net PAR Cab
        Rndir_Car(:,:,j)   = fs * Rnsun_Car(j);                   %   net PAR energy units
        Rndir_PAR(:,:,j)   = fs * Rnsun_PAR(j);                   % [13 x 36 x nl]   net PAR energy units
    end
end

%4.4 total diffuse radiation (net) per leaf area (W m-2 leaf)
for j = 1:nl     % 1 top nl is bottom
    % diffuse incident radiation for the present layer 'j' (mW m-2 um-1)
    E_         = .5*(Emin_(j,:) + Emin_(j+1,:)+ Eplu_(j,:)+ Eplu_(j+1,:));
   
    % incident PAR flux, integrated over all wavelengths (moles m-2 s-1)
    Pdif(j)    = .001 * Sint(e2phot(wlPAR'*1E-9,E_(Ipar),constants),wlPAR);  % [nl] , including conversion mW >> W
  
    % net radiation (mW m-2 um-1) and net PAR (moles m-2 s-1 um-1), per wavelength
    Rndif_(j,:)         = E_.*epsc(j,:);                                                    % [nl,nwl]  Net (absorbed) radiation by leaves
    Pndif_(j,:)         = .001 *(e2phot(wlPAR*1E-9, Rndif_(j,Ipar)',constants))';                     % [nl,nwl]  Net (absorbed) as PAR photons
    Rndif_Cab_(j,:)     = (kChlrel(j,:).*Rndif_(j,spectral.IwlP));    % [nl,nwl]  Net (absorbed) as PAR photons by Cab
    Pndif_Cab_(j,:)     = .001 *(e2phot(spectral.wlP*1E-9, (kChlrel(j,:).*Rndif_(j,spectral.IwlP))',constants))';    % [nl,nwl]  Net (absorbed) as PAR photons by Cab
    Rndif_Car_(j,:)     = (kCarrel(j,:).*Rndif_(j,spectral.IwlP));    % [nl,nwl]  Net (absorbed) as PAR photons by Car
    Pndif_Car_(j,:)     = .001 *(e2phot(spectral.wlP*1E-9, (kCarrel(j,:).*Rndif_(j,spectral.IwlP))',constants))';    % [nl,nwl]  Net (absorbed) as PAR photons by Car

    
    Rndif_PAR_(j,:)     = Rndif_(j,Ipar);                                                   % [nl,nwlPAR]  Net (absorbed) as PAR energy 
   
    % net radiation (W m-2) and net PAR (moles m-2 s-1), integrated over all wavelengths
    Rndif(j)            = .001 * Sint(Rndif_(j,:),wl);              % [nl]  Full spectrum net diffuse flux
    Pndif(j)            =        Sint(Pndif_(j,Ipar),wlPAR);        % [nl]  Absorbed PAR
    Pndif_Cab(j)        =        Sint(Pndif_Cab_(j,:),spectral.wlP);    % [nl]  Absorbed PAR by Cab integrated
    Rndif_Cab(j)        = .001 * Sint(Rndif_Cab_(j,:),spectral.wlP);      % [nl]  Absorbed PAR by Cab integrated
    Pndif_Car(j)        =        Sint(Pndif_Car_(j,:),spectral.wlP);    % [nl]  Absorbed PAR by Car integrated
    Rndif_Car(j)        = .001 * Sint(Rndif_Car_(j,:),spectral.wlP);      % [nl]  Absorbed PAR by Car integrated
    Rndif_PAR(j)        = .001 * Sint(Rndif_PAR_(j,Ipar),wlPAR);    % [nl]  Absorbed PAR by Cab integrated
end

% soil layer, direct and diffuse radiation
Rndirsoil   = .001 * Sint(Esun_.*epss,wl);          % [1] Absorbed solar flux by the soil
Rndifsoil   = .001 * Sint(Emin_(nl+1,:).*epss',wl); % [1] Absorbed diffuse downward flux by the soil (W m-2)

% net (n) radiation R and net PAR P per component: sunlit (u), shaded (h) soil(s) and canopy (c),
% [W m-2 leaf or soil surface um-1]
Rnhc        = Rndif;            % [nl] shaded leaves or needles
Pnhc        = Pndif;            % [nl] shaded leaves or needles
Pnhc_Cab    = Pndif_Cab;        % [nl] shaded leaves or needles
Rnhc_Cab    = Rndif_Cab;        % [nl] shaded leaves or needles
Pnhc_Car    = Pndif_Car;        % [nl] shaded leaves or needles
Rnhc_Car    = Rndif_Car;        % [nl] shaded leaves or needles
Rnhc_PAR    = Rndif_PAR;        % [nl] shaded leaves or needles

if ~options.lite
    [Rnuc,Pnuc,Pnuc_Cab,Rnuc_PAR,Rnuc_Cab,Rnuc_Car,Pnuc_Car] = deal(0*Rndir);
    for j = 1:nl
        %Puc(:,:,j)  = Pdir(:,:,j)      + Pdif(j);      % [13,36,nl] Total fluxes on sunlit leaves or needles
        Rnuc(:,:,j) = Rndir(:,:,j)     + Rndif(j);     % [13,36,nl] Total fluxes on sunlit leaves or needles
        Pnuc(:,:,j)  = Pndir(:,:,j)     + Pndif(j);     % [13,36,nl] Total fluxes on sunlit leaves or needles
        Pnuc_Cab(:,:,j)  = Pndir_Cab(:,:,j)  + Pndif_Cab(j);% [13,36,nl] Total fluxes on sunlit leaves or needles
        Pnuc_Car(:,:,j)  = Pndir_Car(:,:,j)  + Pndif_Car(j);% [13,36,nl] Total fluxes on sunlit leaves or needles
        Rnuc_PAR(:,:,j)  = Rndir_PAR(:,:,j)  + Rndif_PAR(j);% [13,36,nl] Total fluxes on sunlit leaves or needles
        Rnuc_Cab(:,:,j)  = Rndir_Cab(:,:,j)  + Rndif_Cab(j);% [13,36,nl] Total fluxes on sunlit leaves or needles   
        Rnuc_Car(:,:,j)  = Rndir_Car(:,:,j)  + Rndif_Car(j);% [13,36,nl] Total fluxes on sunlit leaves or needles
    end
else   
    %    Puc(:,:,j)  = Pdir      + Pdif(j);      % [13,36,nl] Total fluxes on sunlit leaves or needles
    Rnuc = Rndir     + Rndif;           % [nl] Total fluxes on sunlit leaves or needles
    Pnuc  = Pndir     + Pndif;          % [nl] Total fluxes on sunlit leaves or needles
    Pnuc_Cab  = Pndir_Cab  + Pndif_Cab; % [nl] Total fluxes on sunlit leaves or needles
    Rnuc_Cab  = Rndir_Cab  + Rndif_Cab; % [nl] Total fluxes on sunlit leaves or needles
    Pnuc_Car  = Pndir_Car  + Pndif_Car; % [nl] Total fluxes on sunlit leaves or needles
    Rnuc_Car  = Rndir_Car  + Rndif_Car; % [nl] Total fluxes on sunlit leaves or needles
    Rnuc_PAR  = Rndir_PAR  + Rndif_PAR; % [nl] Total fluxes on sunlit leaves or needles       
end
Rnus        = Rndifsoil + Rndirsoil; % [1] sunlit soil 
Rnhs        = Rndifsoil;  % [1] shaded soil

if options.calc_vert_profiles   
    Pnu1d           = Pnuc;      % [nl]      mean net radiation sunlit leaves
    Pnu1d_Cab       = Pnuc_Cab; 
    Pnu1d_Car       = Pnuc_Car;
    if options.lite == 0
        [Pnu1d  ]           = meanleaf(canopy,Pnuc,         'angles');   % [nli,nlo,nl]      mean net radiation sunlit leaves
        [Pnu1d_Cab  ]       = meanleaf(canopy,Pnuc_Cab,     'angles'); 
        [Pnu1d_Car  ]       = meanleaf(canopy,Pnuc_Car,     'angles');  
    end
    profiles.Pn1d       = ((1-Ps(1:nl)).*Pnhc     + Ps(1:nl).*(Pnu1d));  
    profiles.Pn1d_Cab   = ((1-Ps(1:nl)).*Pnhc_Cab + Ps(1:nl).*(Pnu1d_Cab));  
    profiles.Pn1d_Car   = ((1-Ps(1:nl)).*Pnhc_Car + Ps(1:nl).*(Pnu1d_Car)); 
else
    profiles = struct;
end

%% 5 Model output
% up and down and hemispherical out, cumulative over wavelenght
Eout_       = Eplu_(1,:)';
Eouto       = 0.001 * Sint(Eout_(spectral.IwlP),spectral.wlP);  %     [1] hemispherical out, in optical range (W m-2)
Eoutt       = 0.001 * Sint(Eout_(spectral.IwlT),spectral.wlT);  %     [1] hemispherical out, in thermal range (W m-2)
Lot         = 0.001 * Sint(Lo_(spectral.IwlT),spectral.wlT);    %     [1] hemispherical out, in thermal range (W m-2)

% place output in structure rad
gap.k       = k;        % extinction cofficient in the solar direction
gap.K       = K;        % extinction cofficient in the viewing direction
gap.Ps      = Ps;       % gap fraction in the solar direction
gap.Po      = Po;       % gap fraction in the viewing direction
rad.rsd     = rsd;      % TOC directional-hemispherical reflectance
rad.rdd     = rdd;      % TOC hemispherical-hemispherical reflectance
rad.rdo     = rdo;      % TOC hemispherical-directional reflectance 
rad.rso     = rso;      % TOC directional-directional reflectance   
rad.refl    = Refl;     % TOC reflectance
rad.rho_dd  = rho_dd;   % diffuse-diffuse reflectance for the thin layers
rad.tau_dd  = tau_dd;   % diffuse-diffuse transmittance for the thin layers
rad.rho_sd  = rho_sd;   % direct-diffuse reflectance for the thin layers
rad.tau_ss  = tau_ss;   % direct-direct transmittance for the thin layers
rad.tau_sd  = tau_sd;   % direct-diffuse transmittance for the thin layers
rad.R_sd    = R_sd;
rad.R_dd    = R_dd;
rad.vb      = vb;       % directional backscatter coefficient for diffuse incidence
rad.vf      = vf;       % directional forward scatter coefficient for diffuse incidence
rad.sigf    = sigf;     % forward scatter coefficient for specular flux
rad.sigb    = sigb;     % backscatter coefficient for specular flux
rad.Esun_   = Esun_;    % [nwlx1 double]   incident solar spectrum (mW m-2 um-1)
rad.Esky_   = Esky_;    % [nwlx1 double]   incident sky spectrum (mW m-2 um-1)
rad.PAR     = P*1E6;    % [1 double]       incident spectrally integrated PAR (micromoles m-2 s-1)
rad.EPAR    = EPAR;     % [1 double]       incident PAR in energy units (W m-2)
rad.Eplu_   = Eplu_;    % [nlxnwl double]  upward diffuse radiation in the canopy (mW m-2 um-1)
rad.Emin_   = Emin_;    % [nlxnwl double]  downward diffuse radiation in the canopy (mW m-2 um-1)
rad.Emins_  = Emins_;   % [nlxnwl double]  downward diffuse radiation in the canopy due to direct solar rad (mW m-2 um-1)
rad.Emind_  = Emind_;   % [nlxnwl double]  downward diffuse radiation in the canopy due to sky rad (mW m-2 um-1)
rad.Eplus_  = Eplus_;   % [nlxnwl double]  upward diffuse radiation in the canopy due to direct solar rad (mW m-2 um-1)
rad.Eplud_  = Eplud_;   % [nlxnwl double]  upward diffuse radiation in the canopy due to sky rad (mW m-2 um-1)
rad.Lo_     = Lo_;      % [nwlx1 double]   TOC radiance in observation direction (mW m-2 um-1 sr-1)
rad.Eout_   = Eout_;    % [nwlx1 double]   TOC upward radiation (mW m-2 um-1)
rad.Eouto   = Eouto;    % [1 double]        TOC spectrally integrated upward optical ratiation (W m-2)
rad.Eoutt   = Eoutt;    % [1 double]        TOC spectrally integrated upward thermal ratiation (W m-2)
rad.Lot     = Lot;
rad.Rnhs    = Rnhs;     % [1 double]        net radiation (W m-2) of shaded soil
rad.Rnus    = Rnus;     % [1 double]        net radiation (W m-2) of sunlit soil
rad.Rnhc    = Rnhc;     % [60x1 double]     net radiation (W m-2) of shaded leaves
rad.Rnuc    = Rnuc;     % [13x36x60 double] net radiation (W m-2) of sunlit leaves
rad.Pnh     = 1E6*Pnhc;     % [60x1 double]     net PAR (moles m-2 s-1) of shaded leaves
rad.Pnu     = 1E6*Pnuc;     % [13x36x60 double] net PAR (moles m-2 s-1) of sunlit leaves
rad.Pnh_Cab = 1E6*Pnhc_Cab;% [60x1 double]      net PAR absorbed by Cab (moles m-2 s-1) of shaded leaves
rad.Pnu_Cab = 1E6*Pnuc_Cab; % [13x36x60 double] net PAR absorbed by Cab (moles m-2 s-1) of sunlit leaves
rad.Rnh_Cab = Rnhc_Cab; % [60x1 double]    net PAR absorbed by Cab (W m-2) of shaded leaves
rad.Rnu_Cab = Rnuc_Cab; % [13x36x60 double] net PAR absorbed by Cab (W m-2) of sunlit leaves
rad.Pnh_Car = 1E6*Pnhc_Car;% [60x1 double]      net PAR absorbed by Cab (moles m-2 s-1) of shaded leaves
rad.Pnu_Car = 1E6*Pnuc_Car; % [13x36x60 double] net PAR absorbed by Cab (moles m-2 s-1) of sunlit leaves
rad.Rnh_Car = Rnhc_Car; % [60x1 double]    net PAR absorbed by Cab (W m-2) of shaded leaves
rad.Rnu_Car = Rnuc_Car; % [13x36x60 double] net PAR absorbed by Cab (W m-2) of sunlit leaves
rad.Rnh_PAR = Rnhc_PAR;     % [60x1 double]     net PAR absorbed by Cab (W m-2) of shaded leaves
rad.Rnu_PAR = Rnuc_PAR;     % [13x36x60 double] net PAR absorbed (W m-2) of sunlit
rad.Xdd     =   Xdd;
rad.Xsd     = Xsd;
rad.Xss     = Xss;

% Rn = canopy.LAI*(meanleaf(canopy,rad.Rnhc,'layers',(1-Ps(1:nl)))+meanleaf(canopy,rad.Rnuc,'angles_and_layers',Ps(1:nl)))
% %y1 = canopy.Cv*(rad.Eplu_(end,:)-rad.Eplu_(1,:) + rad.Emin_(1,:) - rad.Emin_(end,:));
% y1 = canopy.Cv*(rad.Eplu_(end,:)-rad.Eplu_(1,:) + rad.Emin_(1,:) - rad.Emin_(end-1,:));
% y2 = canopy.Cv*(1-gap.Ps(end))*rad.Esun_';
% 
% %y2 = (canopy.Cs-gap.Ps(end))*rad.Esun_';
% y = y1+ y2;
% Sint(y,spectral.wlS)*1E-3
% iLAI
%% APPENDIX I function volscat

function [chi_s,chi_o,frho,ftau]    =   volscat(tts,tto,psi,ttli)

%Volscatt version 2.
%created by W. Verhoef
%edited by Joris Timmermans to matlab nomenclature.
% date: 11 February 2008
%tts    [1]         Sun            zenith angle in degrees
%tto    [1]         Observation    zenith angle in degrees
%psi    [1]         Difference of  azimuth angle between solar and viewing position
%ttli   [ttli]      leaf inclination array

deg2rad = pi/180;
nli     = length(ttli);

psi_rad         = psi*deg2rad*ones(nli,1);

cos_psi         = cos(psi*deg2rad);                 %   cosine of relative azimuth angle

cos_ttli        = cos(ttli*deg2rad);                %   cosine of normal of upperside of leaf
sin_ttli        = sin(ttli*deg2rad);                %   sine   of normal of upperside of leaf

cos_tts         = cos(tts*deg2rad);                 %   cosine of sun zenith angle
sin_tts         = sin(tts*deg2rad);                 %   sine   of sun zenith angle

cos_tto         = cos(tto*deg2rad);                 %   cosine of observer zenith angle
sin_tto         = sin(tto*deg2rad);                 %   sine   of observer zenith angle

Cs              = cos_ttli*cos_tts;                 %   p305{1}
Ss              = sin_ttli*sin_tts;                 %   p305{1}

Co              = cos_ttli*cos_tto;                 %   p305{1}
So              = sin_ttli*sin_tto;                 %   p305{1}

As              = max([Ss,Cs],[],2);
Ao              = max([So,Co],[],2);

bts             = acos(-Cs./As);                    %   p305{1}
bto             = acos(-Co./Ao);                    %   p305{2}

chi_o           = 2/pi*((bto-pi/2).*Co + sin(bto).*So);
chi_s           = 2/pi*((bts-pi/2).*Cs + sin(bts).*Ss);

delta1          = abs(bts-bto);                     %   p308{1}
delta2          = pi-abs(bts + bto - pi);           %   p308{1}

Tot             = psi_rad + delta1 + delta2;        %   pag 130{1}

bt1             = min([psi_rad,delta1],[],2);
bt3             = max([psi_rad,delta2],[],2);
bt2             = Tot - bt1 - bt3;

T1              = 2.*Cs.*Co + Ss.*So.*cos_psi;
T2              = sin(bt2).*(2*As.*Ao + Ss.*So.*cos(bt1).*cos(bt3));

Jmin            = (   bt2).*T1 - T2;
Jplus           = (pi-bt2).*T1 + T2;

frho            =  Jplus/(2*pi^2);
ftau            = -Jmin /(2*pi^2);

% pag.309 wl-> pag 135{1}
frho            = max([zeros(nli,1),frho],[],2);
ftau            = max([zeros(nli,1),ftau],[],2);
return

%% APPENDIX II function e2phot

function molphotons = e2phot(lambda,E,constants)
%molphotons = e2phot(lambda,E) calculates the number of moles of photons
%corresponding to E Joules of energy of wavelength lambda (m)

e           = ephoton(lambda,constants);
photons     = E./e;
molphotons  = photons./constants.A;
return;

function E = ephoton(lambda,constants)
%E = phot2e(lambda) calculates the energy content (J) of 1 photon of
%wavelength lambda (m)

h       = constants.h;           % [J s]         Planck's constant
c       = constants.c;           % [m s-1]       speed of light
E       = h*c./lambda;           % [J]           energy of 1 photon
return;

%% APPENDIX III function Pso

function pso    =   Psofunction(K,k,LAI,q,dso,xl)
if dso~=0
    alf         =   (dso/q) *2/(k+K);
    pso         =   exp((K+k)*LAI*xl + sqrt(K*k)*LAI/(alf  )*(1-exp(xl*(alf  ))));% [nl+1]  factor for correlation of Ps and Po
else
    pso         =   exp((K+k)*LAI*xl - sqrt(K*k)*LAI*xl);% [nl+1]  factor for correlation of Ps and Po
    
end
return;

function [R_sd,R_dd,Xss,Xsd,Xdd]  = calc_reflectances(tau_ss,tau_sd,tau_dd,rho_dd,rho_sd,rs,nl,nwl)
[R_sd,R_dd] =   deal(zeros(nl+1,nwl));      % surface reflectance
[Xsd,Xdd]   =   deal(zeros(nl,nwl));        % Effective transmittance
Xss         =   zeros(1,nl);                % Effective transmittance
R_sd(nl+1,:)   =   rs;
R_dd(nl+1,:)   =   rs;
for j=nl:-1:1 % from bottom to top. note nl+1 the background. 1 is the top of canopy.
    Xss      = tau_ss(j);
    dnorm       = 1-rho_dd(j,:).*R_dd(j+1,:);
    Xsd(j,:)    = (tau_sd(j,:)+tau_ss(j).*R_sd(j+1,:).*rho_dd(j,:))./dnorm;
    Xdd(j,:)    = tau_dd(j,:)./dnorm;
    R_sd(j,:)   = rho_sd(j,:)+tau_dd(j,:).*(R_sd(j+1,:).*Xss+R_dd(j+1,:).*Xsd(j,:));
    R_dd(j,:)   = rho_dd(j,:)+tau_dd(j,:).*R_dd(j+1,:).*Xdd(j,:);
end
return;

function [Emin_,Eplu_,Es_] = calc_fluxprofile(Esun_,Esky_,rs,Xss,Xsd,Xdd,R_sd,R_dd,nl,nwl)
[Es_,Emin_,Eplu_]           = deal(zeros(nl+1,nwl));       % [nl+1,nwl]     direct, up and down diff. rad.
Es_(1,:)       = Esun_;
Emin_(1,:)     = Esky_;

for j=1:nl % from top to bottom
    Es_(j+1,:)    =   Xss.*Es_(j,:);
    Emin_(j+1,:)  =   Xsd(j,:).*Es_(j,:)+Xdd(j,:).*Emin_(j,:);
    Eplu_(j,:)    =   R_sd(j,:).*Es_(j,:)+R_dd(j,:).*Emin_(j,:);
end
Eplu_(j+1,:)    =   rs'.*(Es_(j+1,:)+Emin_(j+1,:)); 
%CvdT added calculation of Eplu_(j+1,:)
return;


function [Esun_,Esky_] = calcTOCirr(atmo,meteo,rdd,rsd,wl,nwl)
% calculation of incoming light Esun_ and Esky_
% Extract MODTRAN atmosphere parameters at the SCOPE wavelengths

Fd      = zeros(nwl,1);
Ls      = Planck(wl,meteo.Ta+273.15);

if ~isfield(atmo,'Esun_')
    
    t1  = atmo.M(:,1);
    t3  = atmo.M(:,2);
    t4  = atmo.M(:,3);
    t5  = atmo.M(:,4);
    t12 = atmo.M(:,5);
    t16 = atmo.M(:,6);
    
    % radiation fluxes, downward and upward (these all have dimenstion [nwl]
    % first calculate hemispherical reflectances rsd and rdd according to SAIL
    % these are assumed for the reflectance of the surroundings
    % rdo is computed with SAIL as well
    % assume Fd of surroundings = 0 for the momemnt
    % initial guess of temperature of surroundings from Ta;

    Esun_   = max(1E-6,pi*t1.*t4);
    Esky_   = max(1E-6,pi./(1-t3.*rdd).*(t1.*(t5+t12.*rsd)+Fd+(1-rdd).*Ls.*t3+t16));
    
    % fractional contributions of Esun and Esky to total incident radiation in
    % optical and thermal parts of the spectrum
    if meteo.Rin ~= -999
        % fractional contributions of Esun and Esky to total incident radiation in
        % optical and thermal parts of the spectrum
        
        [fEsuno,fEskyo,fEsunt,fEskyt]          = deal(0*Esun_);   %initialization
        
        J_o             = wl<3000;                          %find optical spectrum
        Esunto          = 0.001 * Sint(Esun_(J_o),wl(J_o)); %Calculate optical sun fluxes (by Integration), including conversion mW >> W
        Eskyto          = 0.001 * Sint(Esky_(J_o),wl(J_o)); %Calculate optical sun fluxes (by Integration)
        Etoto           = Esunto + Eskyto;                  %Calculate total fluxes
        fEsuno(J_o)     = Esun_(J_o)/Etoto;                 %fraction of contribution of Sun fluxes to total light
        fEskyo(J_o)     = Esky_(J_o)/Etoto;                 %fraction of contribution of Sky fluxes to total light
        
        J_t             = wl>=3000;                         %find thermal spectrum
        Esuntt          = 0.001 * Sint(Esun_(J_t),wl(J_t)); %Themal solar fluxes
        Eskytt          = 0.001 * Sint(Esky_(J_t),wl(J_t)); %Thermal Sky fluxes
        Etott           = Eskytt + Esuntt;                  %Total
        fEsunt(J_t)     = Esun_(J_t)/Etott;                 %fraction from Esun
        fEskyt(J_t)     = Esky_(J_t)/Etott;                 %fraction from Esky
        
        Esun_(J_o) = fEsuno(J_o)*meteo.Rin;
        Esky_(J_o) = fEskyo(J_o)*meteo.Rin;
        Esun_(J_t) = fEsunt(J_t)*meteo.Rli;
        Esky_(J_t) = fEskyt(J_t)*meteo.Rli;
    end
    
else
    Esun_ = atmo.Esun_;
    Esky_ = atmo.Esky_;
end
return;

// src/RTMs/RTMt_planck.m

function [rad] = RTMt_planck(spectral,rad,soil,leafopt,canopy,gap,Tcu,Tch,Tsu,Tsh)

% function 'RTMt_sb' calculates total outgoing radiation in hemispherical
% direction and total absorbed radiation per leaf and soil component.
% Radiation is integrated over the whole thermal spectrum with
% Stefan-Boltzman's equation. This function is a simplified version of
% 'RTMt_planck', and is less time consuming since it does not do the
% calculation for each wavelength separately.
%
% Authors: Wout Verhoef and Christiaan van der Tol (c.vandertol@utwente.nl)
% date:     5  Nov 2007
% update:   13 Nov 2007
%           16 Nov 2007 CvdT    improved calculation of net radiation
%           27 Mar 2008 JT      added directional calculation of radiation
%           24 Apr 2008 JT      Introduced dx as thickness of layer (see parameters)
%           31 Oct 2008 JT      introduced optional directional calculation
%           31 Oct 2008 JT      changed initialisation of F1 and F2 -> zeros
%           07 Nov 2008 CvdT    changed layout
%           16 Mar 2009 CvdT    removed Tbright calculation
%              Feb 2013 WV      introduces structures for version 1.40
%           04 Dec 2019 CvdT    adapted for SCOPE-lite
%           17 Mar 2020 CvdT    mSCOPE representation
%
% Table of contents of the function
%   0       preparations
%       0.0     globals
%       0.1     initialisations
%       0.2     parameters
%       0.3     geometric factors of Observer
%       0.4     geometric factors associated with extinction and scattering
%       0.5     geometric factors to be used later with rho and tau
%       0.6     fo for all leaf angle/azumith classes
%   1       calculation of upward and downward fluxes
%   2       total net fluxes
%   Appendix A. Stefan-Boltzmann
%
% usage:
% [rad] = RTMt_sb(constants,rad,soil,leafbio,canopy,gap,Tcu,Tch,Tsu,Tsh,obsdir,spectral)
%
% Most input and output are structures. These structures are further
% specified in a readme file. The temperatures Tcu, Tch, Tsu and Tsh are
% variables.
%
% Input:
%   constants   physical constants
%   rad         a large number of radiative fluxes: spectrally distributed
%               and integrated, and canopy radiative transfer coefficients
%   soil        soil properties
%   leafopt     leaf optical properties
%   canopy      canopy properties (such as LAI and height)
%   gap         probabilities of direct light penetration and viewing
%   Tcu         Temperature of sunlit leaves    (oC), [13x36x60]
%   Tch         Temperature of shaded leaves    (oC), [13x36x60]
%   Tsu         Temperature of sunlit soil      (oC), [1]
%   Tsh         Temperature of shaded soil      (oC), [1]
%
% Output:
%   rad         a large number of radiative fluxes: spectrally distributed
%               and integrated, and canopy radiative transfer coefficients.
%               Here, thermal fluxes are added

%% 0.1 parameters
IT          = spectral.IwlT;   %
wlt         = spectral.wlT;
%deg2rad     = constants.deg2rad;

nl          = canopy.nlayers;
lidf        = canopy.lidf;
Ps          = gap.Ps;
%
rho         = leafopt.refl(:, IT)';    % [1]               Leaf/needle reflection
tau         = leafopt.tran(:, IT)';    % [1]               Leaf/needle transmission
rs          = soil.refl(IT);        % [1]               Soil reflectance
epsc        = 1-rho-tau;              % [nwl]               Emissivity vegetation
epss        = 1-rs;                   % [nwl]               Emissivity soil
LAI         = canopy.LAI;
dx          = 1/nl;
iLAI        = LAI*dx;

Xdd         = rad.Xdd(:,IT);
Xsd         = rad.Xsd(:,IT);
Xss         = repmat(rad.Xss,canopy.nlayers,1);
R_dd        = rad.R_dd(:,IT);
R_sd        = rad.R_sd(:,IT);
rho_dd      = rad.rho_dd(:,IT);
tau_dd      = rad.tau_dd(:,IT);

%% 0.2  initialization of output variables
[piLot_,Eoutte_]    = deal(zeros(1,length(IT))); %          [1,nwlt]
[Emin_,Eplu_]       = deal(zeros(nl+1,length(IT)));       % [nl+1,nwlt]

%% 1. calculation of upward and downward fluxes pag 305

for i = 1:length(IT)
    % 1.1 radiance by components
    Hcsu3          = pi*Planck(wlt(i),Tcu+273.15,epsc(i));
    Hcsh           = pi*Planck(wlt(i),Tch+273.15,epsc(i));
    Hssu           = pi*Planck(wlt(i),Tsu+273.15,epss(i));
    Hssh           = pi*Planck(wlt(i),Tsh+273.15,epss(i));
    % 1.2 radiance by leaf layers Hv and by soil Hs (modified by JAK 2015-01)
    if size(Hcsu3,2)>1
        v1 = repmat( 1/size(Hcsu3, 2), 1, size(Hcsu3, 2)); % vector for computing the mean
        Hcsu2 = reshape(Hcsu3, size(Hcsu3, 1), []);   % create a block matrix from the 3D array
        Hcsu = (v1 * reshape(Hcsu2'*lidf, size(Hcsu3, 2), []))'; % compute column means for each level
    else
        Hcsu = Hcsu3;
    end
    Hc          = Hcsu.*Ps(1:nl) + Hcsh.*(1-Ps(1:nl));      % hemispherical emittance by leaf layers
    Hs          = Hssu.*Ps(nl+1) + Hssh.*(1-Ps(nl+1));      % hemispherical emittance by soil surface
    
    % 1.3 Diffuse radiation
    [U,Es_,Emin,Eplu]           = deal(zeros(nl+1,1));       % [nl+1,nwl]     direct, up and down diff. rad.
    
    U(nl+1)               =   Hs;
    Es_(1)              =   0;
    Emin(1)            =   0;
    
    for j=nl:-1:1      % from bottom to top
        Y(j)  =   (rho_dd(j).*U(j+1)+Hc(j)*iLAI)./(1-rho_dd(j).*R_dd(j+1));
        U(j)  =   tau_dd(j)*(R_dd(j+1).*Y(j)+U(j+1))+Hc(j)*iLAI;
    end
    for j=1:nl       % from top to bottom
        Es_(j+1)    =   Xss(j).*Es_(j);
        Emin(j+1)  =   Xsd(j).*Es_(j)+Xdd(j).*Emin(j)+Y(j);
        Eplu(j)    =   R_sd(j).*Es_(j)+R_dd(j).*Emin(j)+U(j);
    end
    Eplu(nl+1)    =   R_sd(nl).*Es_(nl)+R_dd(nl).*Emin(nl)+Hs;
    
    Emin_(:,i)      = Emin;                        %               downwelling diffuse radiance per layer
    Eplu_(:,i)      = Eplu;                        %               upwelling   diffuse radiance
    
    Eoutte_(i)      = Eplu(1);
    
    % 1.4 Directional radiation and brightness temperature
    K           = gap.K;
    vb          = rad.vb(end);
    vf          = rad.vf(end);
    piLov       = iLAI*...
       (K*Hcsh'*(gap.Po(1:nl)-gap.Pso(1:nl))+  ...              % directional   emitted     radation by shaded leaves
       K*Hcsu'*gap.Pso(1:nl)+ ... % compute column means for each level
       (vb*Emin(1:nl) + vf*Eplu(1:nl))'*gap.Po(1:nl));      % directional   scattered   radiation by vegetation for diffuse incidence   
   
    piLos       = (Hssh*(gap.Po(nl+1)-gap.Pso(nl+1))+ Hssu*gap.Pso(nl+1)); % directional   emitted     radiation by soil
        
    
    piLot_(i)   = piLov + piLos;
    
end
Lot_            = piLot_/pi;

%% 3. Write the output to the rad structure
[rad.Lot_,rad.Eoutte_] = deal(zeros(length(spectral.wlS),1));
rad.Lot_(IT)        = Lot_;
rad.Eoutte_(IT)     = Eoutte_;    %               emitted     diffuse radiance at top
rad.Eplut_          = Eplu_;
rad.Emint_          = Emin_;
return

% 1) CvdT, 11 December 2015.
% We subtract Emin(1), because ALL incident (thermal) radiation from Modtran
% has been taken care of in RTMo. Not ideal but otherwise radiation budget will not close!

// src/RTMs/RTMt_sb.m

function [rad] = RTMt_sb(constants,rad,soil,leafbio,canopy,gap,Tcu,Tch,Tsu,Tsh,obsdir,spectral)

% function 'RTMt_sb' calculates total outgoing radiation in hemispherical
% direction and total absorbed radiation per leaf and soil component.
% Radiation is integrated over the whole thermal spectrum with
% Stefan-Boltzman's equation. This function is a simplified version of
% 'RTMt_planck', and is less time consuming since it does not do the
% calculation for each wavelength separately.
%
% Authors: Wout Verhoef and Christiaan van der Tol (c.vandertol@utwente.nl)
% date:     5  Nov 2007
% update:   13 Nov 2007
%           16 Nov 2007 CvdT    improved calculation of net radiation
%           27 Mar 2008 JT      added directional calculation of radiation
%           24 Apr 2008 JT      Introduced dx as thickness of layer (see parameters)
%           31 Oct 2008 JT      introduced optional directional calculation
%           31 Oct 2008 JT      changed initialisation of F1 and F2 -> zeros
%           07 Nov 2008 CvdT    changed layout
%           16 Mar 2009 CvdT    removed Tbright calculation
%              Feb 2013 WV      introduces structures for version 1.40
%           04 Dec 2019 CvdT    adapted for SCOPE-lite
%           17 Mar 2020 CvdT    mSCOPE representation
%
% Table of contents of the function
%   0       preparations
%       0.0     globals
%       0.1     initialisations
%       0.2     parameters
%       0.3     geometric factors of Observer
%       0.4     geometric factors associated with extinction and scattering
%       0.5     geometric factors to be used later with rho and tau
%       0.6     fo for all leaf angle/azumith classes
%   1       calculation of upward and downward fluxes
%   2       total net fluxes
%   Appendix A. Stefan-Boltzmann
%
% usage:
% [rad] = RTMt_sb(constants,rad,soil,leafbio,canopy,gap,Tcu,Tch,Tsu,Tsh,obsdir,spectral)
%
% Most input and output are structures. These structures are further
% specified in a readme file. The temperatures Tcu, Tch, Tsu and Tsh are
% variables.
%
% Input:
%   constants   physical constants
%   rad         a large number of radiative fluxes: spectrally distributed
%               and integrated, and canopy radiative transfer coefficients
%   soil        soil properties
%   leafopt     leaf optical properties
%   canopy      canopy properties (such as LAI and height)
%   gap         probabilities of direct light penetration and viewing
%   Tcu         Temperature of sunlit leaves    (oC), [13x36x60]
%   Tch         Temperature of shaded leaves    (oC), [13x36x60]
%   Tsu         Temperature of sunlit soil      (oC), [1]
%   Tsh         Temperature of shaded soil      (oC), [1]
%
% Output:
%   rad         a large number of radiative fluxes: spectrally distributed
%               and integrated, and canopy radiative transfer coefficients.
%               Here, thermal fluxes are added

%% 0.1 parameters

nl          = canopy.nlayers;
lidf        = canopy.lidf;
Ps          = gap.Ps;
%
rho         = leafbio.rho_thermal;    % [1]               Leaf/needle reflection
tau         = leafbio.tau_thermal;    % [1]               Leaf/needle transmission
rs          = soil.rs_thermal;        % [1]               Soil reflectance
epsc        = 1-rho-tau;              % [nwl]               Emissivity vegetation
epss        = 1-rs;                   % [nwl]               Emissivity soil
LAI         = canopy.LAI;
dx          = 1/nl;
iLAI        = LAI*dx;

Xdd         = rad.Xdd(:,end);
Xsd         = rad.Xsd(:,end);
Xss         = repmat(rad.Xss,canopy.nlayers,1);
R_dd        = rad.R_dd(:,end);
R_sd        = rad.R_sd(:,end);
rho_dd      = rad.rho_dd(:,end);
tau_dd      = rad.tau_dd(:,end);

%% 1. calculation of upward and downward fluxes pag 305

%1.1 radiance by components
Hcsu3       = epsc*Stefan_Boltzmann(Tcu,constants);%                   Radiance by sunlit leaves
Hcsh        = epsc*Stefan_Boltzmann(Tch,constants);%                   Radiance by shaded leaves
Hssu        = epss*Stefan_Boltzmann(Tsu,constants);%                   Radiance by sunlit soil
Hssh        = epss*Stefan_Boltzmann(Tsh,constants);%                   Radiance by shaded soil

% 1.2 radiance by leaf layers Hv and by soil Hs (modified by JAK 2015-01)
if size(Hcsu3,2)>1
    v1 = repmat( 1/size(Hcsu3, 2), 1, size(Hcsu3, 2)); % vector for computing the mean
    Hcsu2 = reshape(Hcsu3, size(Hcsu3, 1), []);   % create a block matrix from the 3D array
    Hcsu = (v1 * reshape(Hcsu2'*lidf, size(Hcsu3, 2), []))'; % compute column means for each level
else
    Hcsu = Hcsu3;
end
Hc          = Hcsu.*Ps(1:nl) + Hcsh.*(1-Ps(1:nl));      % hemispherical emittance by leaf layers
Hs          = Hssu.*Ps(nl+1) + Hssh.*(1-Ps(nl+1));      % hemispherical emittance by soil surface

% 1.3 Diffuse radiation
[U,Es_,Emin,Eplu]           = deal(zeros(nl+1,1));       % [nl+1,nwl]     direct, up and down diff. rad.

U(nl+1)               =   Hs;
Es_(1)              =   0;
Emin(1)            =   0;

for j=nl:-1:1      % from bottom to top
    Y(j)  =   (rho_dd(j).*U(j+1)+Hc(j)*iLAI)./(1-rho_dd(j).*R_dd(j+1));
    U(j)  =   tau_dd(j)*(R_dd(j+1).*Y(j)+U(j+1))+Hc(j)*iLAI;
end
for j=1:nl       % from top to bottom
    Es_(j+1)    = Xss(j).*Es_(j);
    Emin(j+1)   = Xsd(j).*Es_(j)+Xdd(j).*Emin(j)+Y(j);
    Eplu(j)     = R_sd(j).*Es_(j)+R_dd(j).*Emin(j)+U(j);
end
Eplu(nl+1)      = R_sd(nl).*Es_(nl)+R_dd(nl).*Emin(nl)+Hs;
Eoutte          = Eplu(1);

% 1.4 Directional radiation and brightness temperature
if obsdir
    K           = gap.K;
    vb          = rad.vb(end);
    vf          = rad.vf(end);
    piLov       = iLAI*...
        (K*Hcsh'*(gap.Po(1:nl)-gap.Pso(1:nl))+  ...              % directional   emitted     radation by shaded leaves
        K*Hcsu'*gap.Pso(1:nl)+ ... % compute column means for each level
        (vb*Emin(1:nl) + vf*Eplu(1:nl))'*gap.Po(1:nl));      % directional   scattered   radiation by vegetation for diffuse incidence
    
    piLos       = (Hssh*(gap.Po(nl+1)-gap.Pso(nl+1))+ Hssu*gap.Pso(nl+1));                        % directional   emitted     radiation by  soil

    
    piLot       = piLov + piLos;
    Tbr         = (piLot/constants.sigmaSB)^0.25;
    rad.Lote    = piLot/pi;
    rad.Lot_    = Planck(spectral.wlS,Tbr);% Note that this is the directional blackbody radiance!
    Tbr2        = (Eoutte/constants.sigmaSB)^0.25;
    rad.Eoutte_ = Planck(spectral.wlS,Tbr2);
end

%% 2. total net fluxes
% net radiation per component, in W m-2 (leaf or soil surface)

if size(Hcsu3,2)>1
    Rnuc = 0*Hcsu3;
    for j = 1:nl
        Rnuc(:,:,j) = epsc*(Emin(j) + Eplu(j+1)) - 2*Hcsu3(:,:,j);    % sunlit leaf
    end
else
    Rnuc            = epsc*(Emin(1:end-1) + Eplu(2:end)) - 2*(Hcsu);
end

Rnhc            = epsc*(Emin(1:end-1) + Eplu(2:end)) - 2*(Hcsh);
Rnus            = epss*(Emin(nl+1) - Hssu);                       % sunlit soil
Rnhs            = epss*(Emin(nl+1) - Hssh);                      % shaded soil

%% 3. Write the output to the rad structure
rad.Emint   = Emin;
rad.Eplut   = Eplu;
rad.Eoutte  = Eoutte;
rad.Rnuct   = Rnuc;
rad.Rnhct   = Rnhc;
rad.Rnust   = Rnus;
rad.Rnhst   = Rnhs;
return

% 1) CvdT, 11 December 2015.
% We subtract Emin(1), because ALL incident (thermal) radiation from Modtran
% has been taken care of in RTMo. Not ideal but otherwise radiation budget will not close!

%% Appendix A. Stefan-Boltzmann
function H      =   Stefan_Boltzmann(T_C,constants)

C2K     = constants.C2K;
sigmaSB = constants.sigmaSB;

H       = sigmaSB*(T_C + C2K).^4;
return

// src/RTMs/RTMz.m

function [rad] = RTMz(constants,spectral,rad,soil,leafopt,canopy,gap,angles,Knu,Knh)

% function 'RTMz' calculates the small modification of TOC outgoing
% radiance due to the conversion of Violaxanthin into Zeaxanthin in leaves
%
% Author:  Christiaan van der Tol (c.vandertol@utwente.nl)
% Date:     08 Dec 2016
%           17 Mar 2020     CvdT    added cluming, mSCOPE representation
%           25 Jun 2020     PY      Po, Ps, Pso. fix the problem we have with the oblique angles above 80 degrees

% The inputs and outputs are structures. These structures are further
% specified in a readme file.
%
% Input:
%   spectral    information about wavelengths and resolutions
%   rad         a large number of radiative fluxes: spectrally distributed
%               and integrated, and canopy radiative transfer coefficients.
%   soil        soil properties
%   leafopt     leaf optical properties
%   canopy      canopy properties (such as LAI and height)
%   gap         probabilities of direct light penetration and viewing
%   angles      viewing and observation angles
%
% Output:
%   rad         a large number of radiative fluxes: spectrally distributed
%               and integrated, and canopy radiative transfer coefficients.
%               Here, fluxes are added

%% 0.1 initialisations
wlS         = spectral.wlS';       % SCOPE wavelengths, make column vectors
wlZ         = spectral.wlZ';       % Excitation wavelengths
[dummy,iwlfi]    = intersect(wlS,wlZ); %#ok<ASGLU>
nwlZ        = length(spectral.wlZ);
nl          = canopy.nlayers;
LAI         = canopy.LAI;
iLAI        = LAI/nl;                       % LAI of a layer        [1]
litab       = canopy.litab;
lazitab     = canopy.lazitab;
lidf        = canopy.lidf;
nlazi       = length(lazitab);         % azumith angle
nlinc       = length(litab);           % inclination
nlori       = nlinc * nlazi;           % total number of leaf orientations
layers      = 1:nl;

RZ          = (leafopt.reflZ(:, iwlfi)-leafopt.refl(:, iwlfi))';
TZ          = (leafopt.tranZ(:, iwlfi)-leafopt.tran(:, iwlfi))';

Ps          = gap.Ps;
Po          = gap.Po;
Pso         = gap.Pso;
Qs          = Ps(1:end-1);

% for speed-up the calculation only uses wavelength i and wavelength o part of the spectrum
Esunf_             = rad.Esun_(iwlfi);
[Eminf_,Epluf_]    = deal(zeros(nwlZ,nl+1,2));
Eminf_(:,:,1)      = rad.Emins_(:,iwlfi)';
Eminf_(:,:,2)      = rad.Emind_(:,iwlfi)';
Epluf_(:,:,1)      = rad.Eplus_(:,iwlfi)';
Epluf_(:,:,2)      = rad.Eplud_(:,iwlfi)';

Xdd         = rad.Xdd(:,iwlfi);
rho_dd      = rad.rho_dd(:,iwlfi);
R_dd        = rad.R_dd(:,iwlfi);
tau_dd      = rad.tau_dd(:,iwlfi);
vb          = rad.vb(:,iwlfi);
vf          = rad.vf(:,iwlfi);

%% 0.2 geometric quantities

% geometric factors
deg2rad             = constants.deg2rad;
tto                 = angles.tto;
tts                 = angles.tts;
psi                 = angles.psi;
rs                  = soil.refl(iwlfi,:);           % [nwlfo]     soil reflectance
cos_tto             = cos(tto*deg2rad);             % cos observation zenith angle
sin_tto             = sin(tto*deg2rad);             % sin observation zenith angle

cos_tts             = cos(tts*deg2rad);             % cos solar angle
sin_tts             = sin(tts*deg2rad);             % sin solar angle

cos_ttli            = cos(litab*deg2rad);           % cos leaf inclinaation angles
sin_ttli            = sin(litab*deg2rad);           % sin leaf inclinaation angles
cos_phils           = cos(lazitab*deg2rad);         % cos leaf azimuth angles rel. to sun azi
cos_philo           = cos((lazitab-psi)*deg2rad);   % cos leaf azimuth angles rel. to viewing azi

%% 0.3 geometric factors for all leaf angle/azumith classes
cds                 = cos_ttli*cos_tts*ones(1,36) + sin_ttli*sin_tts*cos_phils;  % [nli,nlazi]
cdo                 = cos_ttli*cos_tto*ones(1,36) + sin_ttli*sin_tto*cos_philo;  % [nli,nlazi]
fs                  = cds/cos_tts;                                               % [nli,nlazi]
absfs               = abs(fs);                                                   % [nli,nlazi]
fo                  = cdo/cos_tto;                                               % [nli,nlazi]
absfo               = abs(fo);                                                   % [nli,nlazi]
fsfo                = fs.*fo;                                                    % [nli,nlazi]
absfsfo             = abs(fsfo);                                                 % [nli,nlazi]
foctl               = fo.*(cos_ttli*ones(1,36));                                 % [nli,nlazi]
fsctl               = fs.*(cos_ttli*ones(1,36));                                 % [nli,nlazi]
ctl2                = cos_ttli.^2*ones(1,36);                                    % [nli,nlazi]

% reshape all the variables
absfs               = reshape(absfs,nlori,1);                                    % [nlori,1]
absfo               = reshape(absfo,nlori,1);                                    % [nlori,1]
fsfo                = reshape(fsfo,nlori,1);                                     % [nlori,1]
absfsfo             = reshape(absfsfo,nlori,1);                                  % [nlori,1]
foctl               = reshape(foctl,nlori,1);                                    % [nlori,1]
fsctl               = reshape(fsctl,nlori,1);                                    % [nlori,1]
ctl2                = reshape(ctl2,nlori,1);                                     % [nlori,1]

%% 1.0 calculation of 'flux' in observation direction
[Fmin_,Fplu_]       = deal(zeros(nl+1,nwlZ,2));
LoF_                = zeros(nwlZ,2);
laz= 1/36;
etah = Kn2Cx(Knh);

if size(Knu,2)==1
    etau = repmat(Kn2Cx(Knu),1,13,36);
else
    etau           = permute(Kn2Cx(Knu),[3 1 2]);    % make dimensions [nl,nlinc,nlazi]
    etau           = reshape(etau,nl,468);    % expand orientations in a vector >> [nl,nlori]
end
etau_lidf = bsxfun(@times,reshape(etau,nlori,nl),repmat(lidf*laz,36,1));     %[nlori,nl]
etah_lidf = bsxfun(@times,repmat(etah,1,nlori)',repmat(lidf*laz,36,1));

for k = 1:2
    [U,Y]     = deal(zeros(nl+1,nwlZ)); 
    MpluEsun  = RZ .* Esunf_*(k<2);      %
    MminEsun  = TZ .* Esunf_*(k<2);
    
    MpluEmin  = RZ  .* Eminf_(:,1:nl,k);	    % [nf,nl+1]  = (nf,ne) * (ne,nl+1)
    MpluEplu  = RZ  .* Epluf_(:,1:nl,k);
    MminEmin  = TZ  .* Eminf_(:,1:nl,k);
    MminEplu  = TZ  .* Epluf_(:,1:nl,k);
    
    wfEs      =  bsxfun(@times,sum(bsxfun(@times,etau_lidf,absfsfo)),MpluEsun) +...
        bsxfun(@times,sum(bsxfun(@times,etau_lidf,fsfo)),MminEsun);
    sfEs     =  bsxfun(@times,sum(bsxfun(@times,etau_lidf,absfs)),MpluEsun) -...
        bsxfun(@times,sum(bsxfun(@times,etau_lidf,fsctl)),MminEsun);
    sbEs     =  bsxfun(@times,sum(bsxfun(@times,etau_lidf,absfs)),MpluEsun) +...
        bsxfun(@times,sum(bsxfun(@times,etau_lidf,fsctl)),MminEsun);
    vfEplu_h  = bsxfun(@times,sum(bsxfun(@times,etah_lidf,absfo)),MpluEplu) -...
        bsxfun(@times,sum(bsxfun(@times,etah_lidf,foctl)),MminEplu);
    vfEplu_u  = bsxfun(@times,sum(bsxfun(@times,etau_lidf,absfo)),MpluEplu) -...
        bsxfun(@times,sum(bsxfun(@times,etau_lidf,foctl)),MminEplu);
    vbEmin_h  = bsxfun(@times,sum(bsxfun(@times,etah_lidf,absfo)),MpluEmin) +...
        bsxfun(@times,sum(bsxfun(@times,etah_lidf,foctl)),MminEmin);
    vbEmin_u  = bsxfun(@times,sum(bsxfun(@times,etau_lidf,absfo)),MpluEmin) +...
        bsxfun(@times,sum(bsxfun(@times,etau_lidf,foctl)),MminEmin);
    sigfEmin_h  = bsxfun(@times,sum(etah_lidf),MpluEmin) -...
        bsxfun(@times,sum(bsxfun(@times,etah_lidf,ctl2)),MminEmin);
    sigfEmin_u  = bsxfun(@times,sum(etau_lidf),MpluEmin) -...
        bsxfun(@times,sum(bsxfun(@times,etau_lidf,ctl2)),MminEmin);
    sigbEmin_h  = bsxfun(@times,sum(etah_lidf),MpluEmin) +...
        bsxfun(@times,sum(bsxfun(@times,etah_lidf,ctl2)),MminEmin);
    sigbEmin_u  = bsxfun(@times,sum(etau_lidf),MpluEmin) +...
        bsxfun(@times,sum(bsxfun(@times,etau_lidf,ctl2)),MminEmin);
    sigfEplu_h  = bsxfun(@times,sum(etah_lidf),MpluEplu) -...
        bsxfun(@times,sum(bsxfun(@times,etah_lidf,ctl2)),MminEplu);
    sigfEplu_u  = bsxfun(@times,sum(etau_lidf),MpluEplu) -...
        bsxfun(@times,sum(bsxfun(@times,etau_lidf,ctl2)),MminEplu);
    sigbEplu_h  = bsxfun(@times,sum(etah_lidf),MpluEplu) +...
        bsxfun(@times,sum(bsxfun(@times,etah_lidf,ctl2)),MminEplu);
    sigbEplu_u  = bsxfun(@times,sum(etau_lidf),MpluEplu) +...
        bsxfun(@times,sum(bsxfun(@times,etau_lidf,ctl2)),MminEplu);
    
    %   Emitted 'flux'
    
    piLs        =   wfEs+vfEplu_u+vbEmin_u;         % sunlit for each layer
    piLd        =   vbEmin_h+vfEplu_h;              % shade leaf for each layer
    Fsmin       =   sfEs+sigfEmin_u+sigbEplu_u;     % Eq. 29a for sunlit leaf
    Fsplu       =   sbEs+sigbEmin_u+sigfEplu_u;     % Eq. 29b for sunlit leaf
    Fdmin       =   sigfEmin_h+sigbEplu_h;          % Eq. 29a for shade leaf
    Fdplu       =   sigbEmin_h+sigfEplu_h;          % Eq. 29b for shade leaf
    Femmin      =   iLAI*bsxfun(@times,Qs', Fsmin) +iLAI* bsxfun(@times,(1-Qs)',Fdmin);
    Femplu      =   iLAI*bsxfun(@times,Qs', Fsplu) +iLAI*bsxfun(@times,(1-Qs)',Fdplu);
    
    for j=nl:-1:1      % from bottom to top
        Y(j,:)  =(rho_dd(j,:).*U(j+1,:)+Femmin(:,j)')./(1-rho_dd(j,:).*R_dd(j+1,:));
        U(j,:) =tau_dd(j,:).*(R_dd(j+1,:).*Y(j,:)+U(j+1,:))+Femplu(:,j)';
    end
    for j=1:nl          % from top to bottom
        Fmin_(j+1,:,k)  = Xdd(j,:).*Fmin_(j,:,k)+Y(j,:);
        Fplu_(j,:,k)  = R_dd(j,:).*Fmin_(j,:,k)+U(j,:);
    end
    piLo1     = iLAI*piLs*Pso(1:nl);
    piLo2     = iLAI*piLd*(Po(1:nl)-Pso(1:nl));
    piLo3     = iLAI*(vb.*Fmin_(layers,:,k)  + vf.*Fplu_(layers,:,k))'*Po(1:nl);
    piLo4     = rs .* (Po(end)* Fmin_(end,:,k)');
    piLtot    = piLo1 + piLo2 + piLo3 + piLo4;
    LoF_(:,k)      = piLtot/pi;
end
Fhem_     = sum(Fplu_(1,:,:),3)';

%% output
rad.Lo_(iwlfi)          = rad.Lo_(iwlfi) + sum(LoF_,2);
rad.rso(iwlfi)          = rad.rso(iwlfi) + LoF_(:,1)./(rad.Esun_(iwlfi));
rad.rdo(iwlfi)          = rad.rdo(iwlfi) + LoF_(:,2)./(rad.Esky_(iwlfi));
rad.refl(iwlfi)         = pi*rad.Lo_(iwlfi)./(rad.Esky_(iwlfi)+rad.Esun_(iwlfi));     % [nwl] 
rad.Eout_(iwlfi)        = rad.Eout_(iwlfi) + Fhem_;

function Cx = Kn2Cx(Kn)
%Cx = 0.70*Kn;  % empirical fit by N Vilfan
Cx = 0.3187*Kn;  % empirical fit by N Vilfan (Vilfan et al, 2018, 2019)
return

// src/RTMs/fluspect_B_CX.m

 function leafopt = fluspect_B_CX(spectral,leafbio,optipar)
%
% function [leafopt] = fluspect(spectral,leafbio,optipar)
% calculates reflectance and transmittance spectra of a leaf using FLUSPECT, 
% plus four excitation-fluorescence matrices
%
% Authors: Wout Verhoef, Christiaan van der Tol (c.vandertol@utwente.nl),
% Joris Timmermans, Nastassia Vilfan
% Date: 2007-2020
% Update from PROSPECT to FLUSPECT: January 2011 (CvdT)
%
%      Nov 2012 (CvdT) Output EF-matrices separately for PSI and PSII
%   31 Jan 2013 (WV)   Adapt to SCOPE v_1.40, using structures for I/O
%   30 May 2013 (WV)   Repair bug in s for non-conservative scattering
%   24 Nov 2013 (WV)   Simplified doubling routine
%   25 Nov 2013 (WV)   Restored piece of code that takes final refl and
%                      tran outputs as a basis for the doubling routine
%   03 Dec 2013 (WV)   Major upgrade. Border interfaces are removed before 
%                      the fluorescence calculation and later added again
%   23 Dec 2013 (WV)   Correct a problem with N = 1 when calculating k 
%                      and s; a test on a = Inf was included
%   01 Apr 2014 (WV)   Add carotenoid concentration (Cca and Kca)
%   19 Jan 2015 (WV)   First beta version for simulation of PRI effect
%   20 Jan 2021 (CvdT) Include PROSPECT-PRO coefficients
% usage:
% [leafopt] = fluspect_b(spectral,leafbio,optipar)
% 
% inputs:
% Cab         = leafbio.Cab;
% Cca         = leafbio.Cca;
% V2Z         = leafbio.V2Z;  % Violaxanthin - Zeaxanthin transition status
%                               [0-1]
% Cw          = leafbio.Cw;
% Cdm         = leafbio.Cdm;
% Cs          = leafbio.Cs;
% N           = leafbio.N; 
% fqe         = leafbio.fqe;
% 
% nr          = optipar.nr;
% Kdm         = optipar.Kdm;
% Kab         = optipar.Kab;
% Kca         = optipar.Kca;
% KcaV        = optipar.KcaV;
% KcaZ        = optipar.KcaZ;
% Kw          = optipar.Kw;
% Ks          = optipar.Ks;
% Kp          = optipar.Kp;
% Kcbc        = optipar.Kcbc;
% phi         = optipar.phi;
% outputs:
% refl          reflectance
% tran          transmittance
% Mb            backward scattering fluorescence matrix, I for PSI and II for PSII
% Mf            forward scattering fluorescence matrix,  I for PSI and II for PSII

%% parameters
% fixed parameters for the fluorescence module
ndub        = 15;           % number of doublings applied
int         = 5;

% Fluspect parameters
Cab         = leafbio.Cab;
Cca         = leafbio.Cca;
V2Z         = leafbio.V2Z;
Cw          = leafbio.Cw;
Cdm         = leafbio.Cdm;
Cs          = leafbio.Cs;
Cant 	    = leafbio.Cant;
Cbc         = leafbio.Cbc;
Cp          = leafbio.Cp;
N           = leafbio.N;
fqe         = leafbio.fqe;

nr          = optipar.nr;
Kdm         = optipar.Kdm;
Kab         = optipar.Kab;

if V2Z == -999 
    % Use old Kca spectrum if this is given as input
    Kca     = optipar.Kca;
else
    % Otherwise make linear combination based on V2Z
    % For V2Z going from 0 to 1 we go from Viola to Zea
    Kca     = (1-V2Z) * optipar.KcaV + V2Z * optipar.KcaZ;    
end
Kw          = optipar.Kw;
Ks          = optipar.Ks;
Kant        = optipar.Kant;
if isfield(optipar, 'Kp')
    Kp        = optipar.Kp;
else
    Kp = 0*Kab;
end
if isfield(optipar, 'Kcbc')
    Kcbc        = optipar.Kcbc;
else
    Kcbc = 0*Kab;
end

phi         = optipar.phi;

%% PROSPECT calculations
Kall        = (Cab*Kab + Cca*Kca + Cdm*Kdm + Cw*Kw  + Cs*Ks + Cant*Kant + Cp*Kp + Cbc*Kcbc)/N;   % Compact leaf layer

j           = find(Kall>0);               % Non-conservative scattering (normal case)
t1          = (1-Kall).*exp(-Kall);
t2          = Kall.^2.*expint(Kall);
tau         = ones(size(t1));
tau(j)      = t1(j)+t2(j);
[kChlrel, kCarrel]     = deal(zeros(size(t1)));
kChlrel(j)  = Cab*Kab(j)./(Kall(j)*N);
kCarrel(j)  = Cca*Kca(j)./(Kall(j)*N);

talf        = calctav(59,nr);
ralf        = 1-talf;
t12         = calctav(90,nr);
r12         = 1-t12;
t21         = t12./(nr.^2);
r21         = 1-t21;

% top surface side
denom       = 1-r21.*r21.*tau.^2;
Ta          = talf.*tau.*t21./denom;
Ra          = ralf+r21.*tau.*Ta;

% bottom surface side
t           = t12.*tau.*t21./denom;
r           = r12+r21.*tau.*t;

% Stokes equations to compute properties of next N-1 layers (N real)
% Normal case

D           = sqrt((1+r+t).*(1+r-t).*(1-r+t).*(1-r-t));
rq          = r.^2;
tq          = t.^2;
a           = (1+rq-tq+D)./(2*r);
b           = (1-rq+tq+D)./(2*t);

bNm1        = b.^(N-1);                  %
bN2         = bNm1.^2;
a2          = a.^2;
denom       = a2.*bN2-1;
Rsub        = a.*(bN2-1)./denom;
Tsub        = bNm1.*(a2-1)./denom;

%			Case of zero absorption
j           = find(r+t >= 1);
Tsub(j)     = t(j)./(t(j)+(1-t(j))*(N-1));
Rsub(j)	    = 1-Tsub(j);

% Reflectance and transmittance of the leaf: combine top layer with next N-1 layers
denom       = 1-Rsub.*r;
tran        = Ta.*Tsub./denom;
refl        = Ra+Ta.*Rsub.*t./denom;

leafopt.refl = refl;
leafopt.tran = tran;
leafopt.kChlrel = kChlrel;
leafopt.kCarrel = kCarrel;
% From here a new path is taken: The doubling method used to calculate
% fluoresence is now only applied to the part of the leaf where absorption
% takes place, that is, the part exclusive of the leaf-air interfaces. The
% reflectance (rho) and transmittance (tau) of this part of the leaf are
% now determined by "subtracting" the interfaces

Rb  = (refl-ralf)./(talf.*t21+(refl-ralf).*r21);  % Remove the top interface
Z   = tran.*(1-Rb.*r21)./(talf.*t21);             % Derive Z from the transmittance

rho = (Rb-r21.*Z.^2)./(1-(r21.*Z).^2);    % Reflectance and transmittance 
tau = (1-Rb.*r21)./(1-(r21.*Z).^2).*Z;    % of the leaf mesophyll layer
t   =   tau;
r   =   max(rho,0);                       % Avoid negative r

% Derive Kubelka-Munk s and k

I_rt     =   (r+t)<1;
D(I_rt)  =   sqrt((1 + r(I_rt) + t(I_rt)) .* ...
                  (1 + r(I_rt) - t(I_rt)) .* ...
                  (1 - r(I_rt) + t(I_rt)) .* ...
                  (1 - r(I_rt) - t(I_rt)));
a(I_rt)  =   (1 + r(I_rt).^2 - t(I_rt).^2 + D(I_rt)) ./ (2*r(I_rt));
b(I_rt)  =   (1 - r(I_rt).^2 + t(I_rt).^2 + D(I_rt)) ./ (2*t(I_rt));    
a(~I_rt) =   1;
b(~I_rt) =   1;

s        =   r./t;
I_a      =   (a>1 & a~=Inf);
s(I_a)   =   2.*a(I_a) ./ (a(I_a).^2 - 1) .* log(b(I_a));

k        =    log(b);
k(I_a)   =   (a(I_a)-1) ./ (a(I_a)+1) .* log(b(I_a));
kChl     =   kChlrel .* k;

%% Fluorescence of the leaf mesophyll layer
% Fluorescence part is skipped for fqe = 0

%light version. The spectral resolution of the irradiance is lowered.
if fqe > 0
    wle         = (400:int:750)';%spectral.wlE';    % excitation wavelengths, transpose to column
    k_iwle      = interp1(spectral.wlP,k,wle);
    s_iwle      = interp1(spectral.wlP,s,wle);
    kChl_iwle   = interp1(spectral.wlP,kChl,wle);
    r21_iwle    = interp1(spectral.wlP,r21,wle);
    rho_iwle    = interp1(spectral.wlP,rho,wle);
    tau_iwle    = interp1(spectral.wlP,tau,wle);
    talf_iwle   = interp1(spectral.wlP,talf,wle);
    
    wlf         = (640:4:850)';%spectral.wlF';    % fluorescence wavelengths, transpose to column
    wlp         = spectral.wlP;     % PROSPECT wavelengths, kept as a row vector

    [~,Iwlf]        = intersect(wlp,wlf);
    eps         = 2^(-ndub);

    % initialisations
    te          = 1-(k_iwle+s_iwle) * eps;   
    tf          = 1-(k(Iwlf)+s(Iwlf)) * eps;  
    re          = s_iwle * eps;
    rf          = s(Iwlf) * eps;

    sigmoid     = 1./(1+exp(-wlf/10)*exp(wle'/10));  % matrix computed as an outproduct
    
    [Mf, Mb] = deal(int*fqe(1) * ((.5*phi(Iwlf))*eps) * kChl_iwle'.*sigmoid);
    
    Ih          = ones(1,length(te));     % row of ones
    Iv          = ones(length(tf),1);     % column of ones

    % Doubling routine
    
    for i = 1:ndub
        
        xe = te./(1-re.*re);  ten = te.*xe;  ren = re.*(1+ten);  
        xf = tf./(1-rf.*rf);  tfn = tf.*xf;  rfn = rf.*(1+tfn);
              
        A11  = xf*Ih + Iv*xe';           A12 = (xf*xe').*(rf*Ih + Iv*re');
        A21  = 1+(xf*xe').*(1+rf*re');   A22 = (xf.*rf)*Ih+Iv*(xe.*re)';
        
        Mfn   = Mf  .* A11 + Mb  .* A12;
        Mbn   = Mb  .* A21 + Mf  .* A22;
        
        te   = ten;  re  = ren;   tf   = tfn;   rf   = rfn;
        Mf  = Mfn; Mb = Mbn;
    end
    
    % Here we add the leaf-air interfaces again for obtaining the final 
    % leaf level fluorescences.
    
    g = Mb; f = Mf;
    
    Rb = rho + tau.^2.*r21./(1-rho.*r21);
    Rb_iwle = interp1(spectral.wlP,Rb,wle);
        
    Xe = Iv * (talf_iwle./(1-r21_iwle.*Rb_iwle))';
    Xf = t21(Iwlf)./(1-r21(Iwlf).*Rb(Iwlf)) * Ih;
    Ye = Iv * (tau_iwle.*r21_iwle./(1-rho_iwle.*r21_iwle))';
    Yf = tau(Iwlf).*r21(Iwlf)./(1-rho(Iwlf).*r21(Iwlf)) * Ih;
    
    A = Xe .* (1 + Ye.*Yf) .* Xf;
    B = Xe .* (Ye + Yf) .* Xf;
    
    gn = A .* g + B .* f;
    fn = A .* f + B .* g;
    
    leafopt.Mb  = gn;
    leafopt.Mf  = fn;
    
end    

return;

function tav = calctav(alfa,nr)

    rd          = pi/180;
    n2          = nr.^2;
    np          = n2+1;
    nm          = n2-1;
    a           = (nr+1).*(nr+1)/2;
    k           = -(n2-1).*(n2-1)/4;
    sa          = sin(alfa.*rd);

    b1          = (alfa~=90)*sqrt((sa.^2-np/2).*(sa.^2-np/2)+k);
    b2          = sa.^2-np/2;
    b           = b1-b2;
    b3          = b.^3;
    a3          = a.^3;
    ts          = (k.^2./(6*b3)+k./b-b/2)-(k.^2./(6*a3)+k./a-a/2);

    tp1         = -2*n2.*(b-a)./(np.^2);
    tp2         = -2*n2.*np.*log(b./a)./(nm.^2);
    tp3         = n2.*(1./b-1./a)/2;
    tp4         = 16*n2.^2.*(n2.^2+1).*log((2*np.*b-nm.^2)./(2*np.*a-nm.^2))./(np.^3.*nm.^2);
    tp5         = 16*n2.^3.*(1./(2*np.*b-nm.^2)-1./(2*np.*a-nm.^2))./(np.^3);
    tp          = tp1+tp2+tp3+tp4+tp5;
    tav         = (ts+tp)./(2*sa.^2);

return;

// src/RTMs/fluspect_mSCOPE.m

function [leafopt]=fluspect_mSCOPE(mly,spectral,leafbio,optipar, nl)
        % leaf reflectance, transmittance and the excitation-fluorescence matrices calculation
        % for 60 sublayers
        indStar =[1,floor(cumsum(mly.pLAI/sum(mly.pLAI))*nl)];  % index of starting for each different layer
        for i=1:mly.nly
            leafbio.Cab     =   mly.pCab(i);
            leafbio.Cw      =   mly.pCw(i);
            leafbio.Cca     =   mly.pCca(i);
            leafbio.Cdm     =   mly.pCdm(i);
            leafbio.Cs      =   mly.pCs(i);
            leafbio.N       =   mly.pN(i);
            [leafopt_ml]    =   fluspect_B_CX(spectral,leafbio,optipar);
            leafopt.refl(i,:)       = leafopt_ml.refl;
            leafopt.tran(i,:)       = leafopt_ml.tran;
            leafopt.Mb(:,:,i)       = leafopt_ml.Mb;
            leafopt.Mf(:,:,i)       = leafopt_ml.Mf;

            leafopt.kChlrel(i,:)    = leafopt_ml.kChlrel;
            leafopt.kCarrel(i,:)    = leafopt_ml.kCarrel;

            in1= indStar(i);
            in2= indStar(i+1);    
            rho_temp(in1:in2,:)    = repmat(leafopt.refl(i,:),in2-in1+1,1);        % [60,nwl]        leaf/needle reflection
            tau_temp(in1:in2,:)    = repmat(leafopt.tran(i,:),in2-in1+1,1);        % [60,nwl]        leaf/needle transmission
            Mb(:,:,in1:in2)        = repmat(leafopt.Mb(:,:,i),[1,1,in2-in1+1]);
            Mf(:,:,in1:in2)        = repmat(leafopt.Mf(:,:,i),[1,1,in2-in1+1]);

            kChlrel_temp(in1:in2,:)= repmat(leafopt.kChlrel(i,:),in2-in1+1,1);
            kCarrel_temp(in1:in2,:)= repmat(leafopt.kCarrel(i,:),in2-in1+1,1);
        end
        leafopt.refl=rho_temp;
        leafopt.tran=tau_temp;
        leafopt.kChlrel = kChlrel_temp;
        leafopt.kCarrel = kCarrel_temp;
        leafopt.Mb     = Mb;
        leafopt.Mf     = Mf;

end

// src/fluxes/biochemical.m

function biochem_out = biochemical(leafbio,meteo,options,constants,fV)
%
% Date: 	21 Sep 2012
% Update:   20 Feb 2013
% Update:      Aug 2013: correction of L171: Ci = Ci*1e6 ./ p .* 1E3;
% Update:   2016-10 - (JAK) major rewrite to accomodate an iterative solution to the Ball-Berry equation
%                   - also allows for g_m to be specified for C3 plants, but only if Ci_input is provided.
% Update:   25 Feb 2021: Temperature reponse functions by Dutta et al.
%                       implemented
% Authors: 	Joe Berry and Christiaan van der Tol, Ari Kornfeld, contributions of others.
% Sources:
%           Farquhar et al. 1980, Collatz et al (1991, 1992), and:
%
% Dutta, D., Schimel, D. S., Sun, Y., Tol, C. V. D., & Frankenberg, C. (2019). 
% Optimal inverse estimation of ecosystem parameters from observations of carbon and energy fluxes. 
% Biogeosciences, 16(1), 77-103.
%
% Van der Tol, C., Berry, J. A., Campbell, P. K. E., & Rascher, U. (2014). 
% Models of fluorescence and photosynthesis for interpreting measurements 
% of solar?induced chlorophyll fluorescence. 
% Journal of Geophysical Research: Biogeosciences, 119(12), 2312-2327.
%
% Bonan, G. B., Lawrence, P. J., Oleson, K. W., Levis, S., Jung, M., 
% Reichstein, M., ... & Swenson, S. C. (2011). 
% Improving canopy processes in the Community Land Model version 4 (CLM4) 
% using global flux fields empirically inferred from FLUXNET data. 
% Journal of Geophysical Research: Biogeosciences, 116(G2).
%
%
% This function calculates:
%    - stomatal resistance of a leaf or needle (s m-1)
%    - photosynthesis of a leaf or needle (umol m-2 s-1)
%    - fluorescence of a leaf or needle (fraction of fluor. in the dark)
%
% Usage:
% biochem_out = biochemical(leafbio,meteo,options,constants,fV)
% the function was tested for Matlab R2017a
%
% Calculates net assimilation rate A, fluorescence F using biochemical model
%
% Input (units are important):
% structure 'leafbio' with the following elements:
% Knparams   % [], [], []           parameters for empirical Kn (NPQ) model: Kn = Kno * (1+beta).*x.^alpha./(beta + x.^alpha);
%       [Kno, Kn_alpha, Kn_beta]
%  or, better, as individual fields:
%   Kno                                     Kno - the maximum Kn value ("high light")
%   Kn_alpha, Kn_beta                      alpha, beta: curvature parameters
%
% Cs        % [ppmV or umol mol]    initial estimate of conc. of CO2 in the
%                                   ...bounary layer of the leaf
% Q         % [umol photons m-2 s-1]net radiation, PAR
% fPAR     % [0-1]                 fraction of incident light that is absorbed by the leaf (default = 1, for compatibility)
% T         % [oC or K]             leaf temperature
% eb        % [hPa = mbar]          intial estimate of the vapour pressure in leaf boundary layer
% O         % [mmol/mol]            concentration of O2 (in the boundary
%                                   ...layer, but no problem to use ambient)
% p         % [hPa]                 air pressure
% Vcmax25 (Vcmo)  % [umol/m2/s]     maximum carboxylation capacity @ 25 degC
% BallBerrySlope (m) % []           Ball-Berry coefficient 'm' for stomatal regulation
% BallBerry0 % []              (OPTIONAL) Ball-Berry intercept term 'b' (if present, an iterative solution is used)
%                                     setting this to zeo disables iteration. Default = 0.01
%
% Type      % ['C3', 'C4']          text parameter, either 'C3' for C3 or any
%                                   ...other text for C4
% tempcor   % [0, 1]               boolean (0 or 1) whether or not
%                                   ...temperature correction to Vcmax has to be applied.
%
% effcon    [mol CO2/mol e-]  number of CO2 per electrons - typically 1/5 for C3 and 1/6 for C4

% RdPerVcmax25 (Rdparam)  % []     respiration as fraction of Vcmax25
% stressfactor [0-1]               stress factor to reduce Vcmax (for
%                                   example soil moisture, leaf age). Use 1 to "disable" (1 = no stress)
%  OPTIONAL
% Kpep25 (kp)   % [umol/m2/s]         PEPcase activity at 25 deg C (defaults to Vcmax/56
% atheta      % [0-1]                  smoothing parameter for transition between Vc and Ve (light- and carboxylation-limited photosynthesis)
% useTLforC3  % boolean              whether to enable low-temperature attenuation of Vcmax in C3 plants (its always on for C4 plants)
% po0         %  double            Kp,0 (Kp,max) = Fv/Fm (for curve fitting)
% g_m         % mol/m2/s/bar      Mesophyll conductance (default: Infinity, i.e. no effect of g_m)

% Note: always use the prescribed units. Temperature can be either oC or K
% Note: input can be single numbers, vectors, or n-dimensional
% matrices
%
% Output:
% structure 'biochem_out' with the following elements:
% A         % [umol/m2/s]           net assimilation rate of the leaves
% Cs        % [umol/m3]             CO2 concentration in the boundary layer
% eta0      % []                    fluorescence as fraction of dark
%                                   ...adapted (fs/fo)
% rcw       % [s m-1]               stomatal resistance
% qE        % []                    non photochemical quenching
% fs        % []                    fluorescence as fraction of PAR
% Ci        % [umol/m3]             internal CO2 concentration
% Kn        % []                    rate constant for excess heat
% fo        % []                    dark adapted fluorescence (fraction of aPAR)
% fm        % []                    light saturated fluorescence (fraction of aPAR)
% qQ        % []                    photochemical quenching
% Vcmax     % [umol/m2/s]           carboxylation capacity after
%                                   ... temperature correction

%% physiological options
tempcor       = options.apply_T_corr;

%% input

% atmospheric physics
rhoa            = constants.rhoa;           % [kg m-3]       specific mass of air
Mair            = constants.Mair;           % [g mol-1]      molecular mass of dry air
R               = constants.R;              % [J mol-1K-1]   Molar gas constant
Q               = meteo.Q;                  % [umol m-2 s-1] absorbed PAR flux
Cs              = meteo.Cs;
T               = meteo.T + 273.15*(meteo.T<200); % convert temperatures to K if not already
eb              = meteo.eb;
O               = meteo.Oa;
p               = meteo.p;

%biochemical
Type            = leafbio.Type;
Vcmax25         = fV.*leafbio.Vcmax25;
BallBerrySlope  = leafbio.BallBerrySlope;
RdPerVcmax25    = leafbio.RdPerVcmax25;
BallBerry0      = leafbio.BallBerry0;
Tref            = 25+273.15;        % [K]           absolute temperature at 25 oC

%Jmax25      = 1.97*Vcmax25;
%TPU25           = 0.06*Jmax25;              % Triose Phosphate utilization rate

% All of the Following Values are Adopted from Bonan et. al 2011 paper and
% pertaings to C3 photosythesis

Kc25            = 405;                   % [umol mol-1]
Ko25            = 279;                   % [mmol mol-1]
spfy25          = 2444;                  % specificity (Computed from Bernacchhi et al 2001 paper)

% convert all to bar: CO2 was supplied in ppm, O2 in permil, and pressure in mBar
ppm2bar     =  1e-6 .* (p .*1E-3);
Cs          = Cs .* ppm2bar;
O           = (O * 1e-3) .* (p .*1E-3) .* strcmp('C3',Type);    % force O to be zero for C4 vegetation (this is a trick to prevent oxygenase)
Kc25        = Kc25 * 1e-6;
Ko25        = Ko25 * 1e-3;
Gamma_star25    = 0.5 .*O./spfy25;      % [ppm] compensation point in absence of Rd
Rd25            = RdPerVcmax25 * Vcmax25;
if strcmpi('C3', Type)
    effcon =  1/5;
else
    effcon = 1/6; % C4
end
atheta      = 0.8;

% Mesophyll conductance: by default we ignore its effect
%  so Cc = Ci - A/gm = Ci
g_m = Inf;
if isfield(leafbio, 'g_m')
    g_m = leafbio.g_m * 1e6; % convert from mol to umol
end
stressfactor  = leafbio.stressfactor;

% fluorescence
Knparams    = [leafbio.Kn0, leafbio.Knalpha, leafbio.Knbeta];
Kf          = 0.05;             % []            rate constant for fluorescence
Kd          = max(0.8738,  0.0301*(T-273.15)+ 0.0773);
Kp          = 4.0;              % []            rate constant for photochemisty

%% temperature corrections

[f.Vcmax,f.Rd,f.TPU,f.Kc,f.Ko,f.Gamma_star] = deal(1);
if tempcor
    if strcmpi('C4', Type)
        % RdPerVcmax25 = 0.025;  % Rd25 for C4 is different than C3
        %   Rd25 = RdPerVcmax25 * Vcmax25;
        % Constant parameters for temperature correction of Vcmax
        Q10 = leafbio.TDP.Q10;                           % Unit is  []
        s1  = leafbio.TDP.s1;                            % Unit is [K]
        s2  = leafbio.TDP.s2;                            % Unit is [K^-1]
        s3  = leafbio.TDP.s3;                            % Unit is [K]
        s4  = leafbio.TDP.s4;                            % Unit is [K^-1]
        
        % Constant parameters for temperature correction of Rd
        s5  = leafbio.TDP.s5;                            % Unit is [K]
        s6  = leafbio.TDP.s6;                            % Unit is [K^-1]
        
        fHTv = 1 + exp(s1.*(T - s2));
        fLTv = 1 + exp(s3.*(s4 - T));
        Vcmax = (Vcmax25 .* Q10.^(0.1.*(T-Tref)))./(fHTv .* fLTv); % Temp Corrected Vcmax
        
        % Temperature correction of Rd
        
        fHTv = 1 + exp(s5.*(T - s6));
        Rd = (Rd25 .* Q10.^(0.1.*(T-Tref)))./fHTv; % Temp Corrected Rd
        % Temperature correction of Ke
        Ke25 = 20000 .* Vcmax25 ;               % Unit is  []
        Ke = (Ke25 .* Q10.^(0.1.*(T-Tref)));    % Temp Corrected Ke
        
    else
        % temperature correction of Vcmax
        deltaHa     = leafbio.TDP.delHaV;                % Unit is  [J K^-1]
        deltaS      = leafbio.TDP.delSV;                 % unit is [J mol^-1 K^-1]
        deltaHd     = leafbio.TDP.delHdV;                % unit is [J mol^-1]
        fTv         = temperature_functionC3(Tref,R,T,deltaHa);
        fHTv        = high_temp_inhibtionC3(Tref,R,T,deltaS,deltaHd);
        f.Vcmax     = fTv .* fHTv;
        
%         % temperature correction for TPU
%         deltaHa     = leafbio.TDP.delHaP;                % Unit is  [J K^-1]
%         deltaS      = leafbio.TDP.delSP;                 % unit is [J mol^-1 K^-1]
%         deltaHd     = leafbio.TDP.delHdP;                % unit is [J mol^-1]
%         fTv         = temperature_functionC3(Tref,R,T,deltaHa);
%         fHTv        = high_temp_inhibtionC3(Tref,R,T,deltaS,deltaHd);
%         f.TPU       = fTv .* fHTv;
        
        % temperature correction for Rd
        deltaHa     = leafbio.TDP.delHaR;               % Unit is  [J K^-1]
        deltaS      = leafbio.TDP.delSR;                % unit is [J mol^-1 K^-1]
        deltaHd     = leafbio.TDP.delHdR;               % unit is [J mol^-1]
        fTv         = temperature_functionC3(Tref,R,T,deltaHa);
        fHTv        = high_temp_inhibtionC3(Tref,R,T,deltaS,deltaHd);
        f.Rd        = fTv .* fHTv;
        
        % temperature correction for Kc
        deltaHa     = leafbio.TDP.delHaKc;               % Unit is  [J K^-1]
        fTv         = temperature_functionC3(Tref,R,T,deltaHa);
        f.Kc        = fTv;
        
        % temperature correction for Ko
        deltaHa     = leafbio.TDP.delHaKo;               % Unit is  [J K^-1]
        fTv         = temperature_functionC3(Tref,R,T,deltaHa);
        f.Ko        = fTv;
        
        % temperature correction for Gamma_star
        deltaHa     = leafbio.TDP.delHaT;               % Unit is  [J K^-1]
        fTv         = temperature_functionC3(Tref,R,T,deltaHa);
        f.Gamma_star = fTv;
        
        Ke          = 1; % dummy value (only needed for C4)
    end
else
    Ke          = 1; % dummy value (only needed for C4)
end

if strcmp('C3', Type)
    Vcmax       = Vcmax25.*f.Vcmax.*stressfactor;
    Rd          = Rd25    .*f.Rd.* stressfactor;
    %TPU         = TPU25   .* f.TPU;
    Kc          = Kc25    .* f.Kc;
    Ko          = Ko25    .* f.Ko;
end
Gamma_star   = Gamma_star25 .* f.Gamma_star;

%% calculation of potential electron transport rate
po0         = Kp./(Kf+Kd+Kp);         % maximum dark photochemistry fraction, i.e. Kn = 0 (Genty et al., 1989)
Je          = 0.5*po0 .* Q;          % potential electron transport rate (JAK: add fPAR);
%Gamma_star  = 0.5 .*O ./spfy; %[bar]       compensation point in absence of Rd (i.e. gamma*) [bar]

if strcmp(Type, 'C3')
    MM_consts = (Kc .* (1+O./Ko)); % Michaelis-Menten constants
    Vs_C3 = (Vcmax/2);
    %  minimum Ci (as fraction of Cs) for BallBerry Ci. (If Ci_input is present we need this only as a placeholder for the function call)
    minCi = 0.3;
else
    % C4
    MM_consts = 0; % just for formality, so MM_consts is initialized
    Vs_C3 = 0;     %  the same
    minCi = 0.1;  % C4
end

%% calculation of Ci (internal CO2 concentration)
RH = min(1, eb./satvap(T-273.15) ); % jak: don't allow "supersaturated" air! (esp. on T curves)
computeA()  % clears persistent fcount
computeA_fun = @(x) computeA(x, Type, g_m, Vs_C3, MM_consts, Rd, Vcmax, Gamma_star, Je, effcon, atheta, Ke);

if all(BallBerry0 == 0)
    % b = 0: no need to iterate:
    Ci = BallBerry(Cs, RH, [], BallBerrySlope, BallBerry0, minCi);
    %     A =  computeA_fun(Ci);   
else
    % compute Ci using iteration (JAK)
    % it would be nice to use a built-in root-seeking function but fzero requires scalar inputs and outputs,
    % Here I use a fully vectorized method based on Brent's method (like fzero) with some optimizations.
    tol = 1e-7;  % 0.1 ppm more-or-less
    % Setting the "corner" argument to Gamma may be useful for low Ci cases, but not very useful for atmospheric CO2, so it's ignored.
    %                     (fn,                           x0, corner, tolerance)
    [Ci] = fixedp_brent_ari(@(x) Ci_next(x, Cs, RH, minCi, BallBerrySlope, BallBerry0, computeA_fun, ppm2bar), Cs, [], tol); % [] in place of Gamma: it didn't make much difference
    %NOTE: A is computed in Ci_next on the final returned Ci. fixedp_brent_ari() guarantees that it was done on the returned values.
    %     A =  computeA_fun(Ci);
end

[A, biochem_out]    = computeA_fun(Ci);
Ag                  = biochem_out.Ag;
CO2_per_electron    = biochem_out.CO2_per_electron;

%% Compute A, electron transport rate, and stomatal resistance
gs                  = max(0, 1.6 * A .* ppm2bar ./ (Cs-Ci));     % stomatal conductance
Ja                  = Ag ./ CO2_per_electron;   % actual electron transport rate
rcw                 = (rhoa./(Mair*1E-3))./gs;    % stomatal resistance

%% fluorescence
ps          = po0.*Ja./Je;               % this is the photochemical yield
nanPs = isnan(ps);
if any(nanPs)
    if numel(po0) == 1
        ps(nanPs) = po0;
    else
        ps(nanPs) = po0(nanPs);  % happens when Q = 0, so ps = po0 (other cases of NaN have been resolved)
    end
end
ps_rel   = max(0,  1-ps./po0);       % degree of light saturation: 'x' (van der Tol e.a. 2014)

[eta,qE,qQ,fs,fo,fm,fo0,fm0,Kn]    = Fluorescencemodel(ps, ps_rel, Kp,Kf,Kd,Knparams);
Kpa         = ps./fs*Kf;

%% convert back to ppm
Cc = [];
if ~isempty(g_m)
    Cc    = (Ci - A/g_m) ./ ppm2bar;
end
Ci          = Ci  ./ ppm2bar;
%Cs          = Cs  ./ ppm2bar;

%% Collect outputs
biochem_out.A       = A;
biochem_out.Ci      = Ci;
if ~isempty(Cc)
    biochem_out.Cc = Cc;
end
biochem_out.rcw     = rcw;
biochem_out.gs      =  gs;
biochem_out.RH      =  RH;
biochem_out.Vcmax   = Vcmax;
biochem_out.Rd      = Rd;
biochem_out.Ja      = Ja;
biochem_out.ps      = ps; % photochemical yield
biochem_out.ps_rel  = ps_rel;   % degree of ETR saturation 'x' (van der Tol e.a. 2014)
biochem_out.Kd      = Kd;  % K_dark(T)
biochem_out.Kn      = Kn;  % K_n(x);  x = 1 - ps/p00 == 1 - Ja/Je
biochem_out.NPQ     = Kn ./ (Kf + Kd); % 
biochem_out.Kf      = Kf;  % Kf = 0.05 (const)
biochem_out.Kp0     = Kp;  % Kp = 4.0 (const): Kp, max
biochem_out.Kp      = Kpa; % Kp,actual
biochem_out.eta     = eta;
biochem_out.qE      = qE;
biochem_out.fs      = fs;  % keep this for compatibility with SCOPE
biochem_out.ft      = fs;  % keep this 
biochem_out.SIF     = fs .* Q;
biochem_out.fo0     = fo0;
biochem_out.fm0     = fm0;
biochem_out.fo      = fo;
biochem_out.fm      = fm;
biochem_out.Fm_Fo   = fm ./ fo;  % parameters used for curve fitting
biochem_out.Ft_Fo   = fs ./ fo;  % parameters used for curve fitting
biochem_out.qQ      = qQ;
biochem_out.Phi_N   = Kn./(Kn +Kp+Kf+Kd);
return;

end  % end of function biochemical


%% quadratic formula, root of least magnitude
function x = sel_root(a,b,c, dsign)
%  sel_root - select a root based on the fourth arg (dsign = discriminant sign)
%    for the eqn ax^2 + bx + c,
%    if dsign is:
%       -1, 0: choose the smaller root
%       +1: choose the larger root
%  NOTE: technically, we should check a, but in biochemical, a is always > 0
if a == 0  % note: this works because 'a' is a scalar parameter!
    x      = -c./b;
else
    if any(dsign == 0)
        dsign(dsign == 0) = -1; % technically, dsign==0 iff b = c = 0, so this isn't strictly necessary except, possibly for ill-formed cases)
    end
    %disc_root = sqrt(b.^2 - 4.*a.*c); % square root of the discriminant (doesn't need a separate line anymore)
    %  in MATLAB (2013b) assigning the intermediate variable actually slows down the code! (~25%)
    x = (-b + dsign.* sqrt(b.^2 - 4.*a.*c))./(2.*a);
end
end %of min_root of quadratic formula


%% Ball Berry Model
function [Ci, gs] = BallBerry(Cs, RH, A, BallBerrySlope, BallBerry0, minCi, Ci_input)
%  Cs  : CO2 at leaf surface
%  RH  : relative humidity
%  A   : Net assimilation in 'same units of CO2 as Cs'/m2/s
% BallBerrySlope, BallBerry0,
% minCi : minimum Ci as a fraction of Cs (in case RH is very low?)
% Ci_input : will calculate gs if A is specified.
if nargin > 6 && ~isempty(Ci_input)
    % Ci is given: try and compute gs
    Ci = Ci_input;
    gs = [];
    if ~isempty(A) && nargout > 1
        gs = gsFun(Cs, RH, A, BallBerrySlope, BallBerry0);
    end
elseif all(BallBerry0 == 0) || isempty(A)
    % EXPLANATION:   *at equilibrium* CO2_in = CO2_out => A = gs(Cs - Ci) [1]
    %  so Ci = Cs - A/gs (at equilibrium)                                 [2]
    %  Ball-Berry suggest: gs = m (A RH)/Cs + b   (also at equilib., see Leuning 1990)
    %  if b = 0 we can rearrange B-B for the second term in [2]:  A/gs = Cs/(m RH)
    %  Substituting into [2]
    %  Ci = Cs - Cs/(m RH) = Cs ( 1- 1/(m RH)  [ the 1.6 converts from CO2- to H2O-diffusion ]
    Ci      = max(minCi .* Cs,  Cs.*(1-1.6./(BallBerrySlope .* RH)));
    gs = [];
else
    %  if b > 0  Ci = Cs( 1 - 1/(m RH + b Cs/A) )
    % if we use Leuning 1990, Ci = Cs - (Cs - Gamma)/(m RH + b(Cs - Gamma)/A)  [see def of Gamma, above]
    % note: the original B-B units are A: umol/m2/s, ci ppm (umol/mol), RH (unitless)
    %   Cs input was ppm but was multiplied by ppm2bar above, so multiply A by ppm2bar to put them on the same scale.
    %  don't let gs go below its minimum value (i.e. when A goes negative)
    gs = gsFun(Cs, RH, A, BallBerrySlope, BallBerry0);
    Ci = max(minCi .* Cs,  Cs - 1.6 * A./gs) ;
end

end % function

function gs = gsFun(Cs, RH, A, BallBerrySlope, BallBerry0)
% add in a bit just to avoid div zero. 1 ppm = 1e-6 (note since A < 0 if Cs ==0, it gives a small gs rather than maximal gs
gs = max(BallBerry0,  BallBerrySlope.* A .* RH ./ (Cs+1e-9)  + BallBerry0);
% clean it up:
%gs( Cs == 0 ) = would need to be max gs here;  % eliminate infinities
gs( isnan(Cs) ) = NaN;  % max(NaN, X) = X  (MATLAB 2013b) so fix it here
end


%% Fluorescence model
function [eta,qE,qQ,fs,fo,fm,fo0,fm0,Kn] = Fluorescencemodel(ps,x, Kp,Kf,Kd,Knparams)
% note: x isn't strictly needed as an input parameter but it avoids code-duplication (of po0) and it's inherent risks.

Kno = Knparams(1);
alpha = Knparams(2);
beta = Knparams(3);

x_alpha = exp(log(x).*alpha); % this is the most expensive operation in this fn; doing it twice almost doubles the time spent here (MATLAB 2013b doesn't optimize the duplicate code)
Kn = Kno * (1+beta).* x_alpha./(beta + x_alpha);

fo0         = Kf./(Kf+Kp+Kd);        % dark-adapted fluorescence yield Fo,0
fo          = Kf./(Kf+Kp+Kd+Kn);     % light-adapted fluorescence yield in the dark Fo
fm          = Kf./(Kf   +Kd+Kn);     % light-adapted fluorescence yield Fm
fm0         = Kf./(Kf   +Kd);        % dark-adapted fluorescence yield Fm
fs          = fm.*(1-ps);            % steady-state (light-adapted) yield Ft (aka Fs)
eta         = fs./fo0;
qQ          = 1-(fs-fo)./(fm-fo);    % photochemical quenching
qE          = 1-(fm-fo)./(fm0-fo0);  % non-photochemical quenching

end

%% Test-function for iteration
%   (note that it assigns A in the function's context.)
%   As with the next section, this code can be read as if the function body executed at this point.
%    (if iteration was used). In other words, A is assigned at this point in the file (when iterating).
function [err, Ci_out] = Ci_next(Ci_in, Cs, RH, minCi, BallBerrySlope, BallBerry0, A_fun, ppm2bar)
% compute the difference between "guessed" Ci (Ci_in) and Ci computed using BB after computing A
A = A_fun(Ci_in);
A_bar = A .* ppm2bar;
Ci_out = BallBerry(Cs, RH, A_bar, BallBerrySlope, BallBerry0, minCi); %[Ci_out, gs]

err = Ci_out - Ci_in; % f(x) - x
end

%% Compute Assimilation.
%  Note: even though computeA() is written as a separate function,
%    the code is, in fact, executed exactly this point in the file (i.e. between the previous if clause and the next section
function [A, biochem_out] = computeA(Ci, Type, g_m, Vs_C3, MM_consts, Rd, Vcmax, Gamma_star, Je, effcon, atheta, kpepcase)
% global: Type, Vcmax, Gamma_star, MM_consts, Vs_C3, effcon, Je, atheta, Rd    %Kc, O, Ko, Vcmax25, qt
persistent fcount
if nargin == 0
    fcount = 0;
    return
end
if strcmpi('C3', Type)
    %[Ci, gs] = BallBerry(Cs, RH, A_bar, BallBerrySlope, BallBerry0, 0.3, Ci_input);
    %effcon      = 0.2;
    % without g_m:
    Vs          = Vs_C3; % = (Vcmax25/2) .* exp(log(1.8).*qt);    % doesn't change on iteration.
    if any(g_m < Inf)
        % with g_m:
        Vc = sel_root( 1./g_m, -(MM_consts + Ci +(Rd + Vcmax)./g_m), Vcmax.*(Ci - Gamma_star + Rd./g_m), -1);
        Ve = sel_root( 1./g_m, -(Ci + 2*Gamma_star +(Rd + Je .* effcon)./g_m), Je .* effcon.*(Ci - Gamma_star + Rd./g_m), -1);
        CO2_per_electron = Ve ./ Je;
    else
        Vc          = Vcmax.*(Ci-Gamma_star)./(MM_consts + Ci);  % MM_consts = (Kc .* (1+O./Ko)) % doesn't change on iteration.
        CO2_per_electron = (Ci-Gamma_star)./(Ci+2*Gamma_star) .* effcon;
        Ve          = Je .* CO2_per_electron;
    end
else  %C4
    %[Ci, gs] = BallBerry(Cs, RH, A_bar, BallBerrySlope, BallBerry0, 0.1, Ci_input);
    Vc          = Vcmax;
    Vs          = kpepcase.*Ci;
    %effcon      = 0.17;                    % Berry and Farquhar (1978): 1/0.167 = 6
    CO2_per_electron = effcon; % note: (Ci-Gamma_star)./(Ci+2*Gamma_star) = 1 for C4 (since O = 0); this line avoids 0/0 when Ci = 0
    Ve          = Je .* CO2_per_electron;
end

% find the smoothed minimum of Ve, Vc = V, then V, Vs
%         [a1,a2]     = abc(atheta,-(Vc+Ve),Vc.*Ve);
%         % select the min or max  depending on the side of the CO2 compensation point
%         %  note that Vc, Ve < 0 when Ci < Gamma_star (as long as Q > 0; Q = 0 is also ok),
%         %     so the original construction selects the value closest to zero.
%         V           = min(a1,a2).*(Ci>Gamma_star) + max(a1,a2).*(Ci<=Gamma_star);
%         [a1,a2]     = abc(0.98,-(V+Vs),V.*Vs);
%         Ag          = min(a1,a2);
V           = sel_root(atheta,-(Vc+Ve),Vc.*Ve, sign(-Vc) ); % i.e. sign(Gamma_star - Ci)
Ag          = sel_root(0.98,-(V+Vs),V.*Vs, -1);
A           = Ag - Rd;
fcount = fcount + 1; % # of times we called computeA

if nargout > 1
    biochem_out.A = A;
    biochem_out.Ag = Ag;
    biochem_out.Vc = Vc;
    biochem_out.Vs = Vs;
    biochem_out.Ve = Ve;
    biochem_out.CO2_per_electron = CO2_per_electron;
    biochem_out.fcount = fcount;
end

end

%% Temperature Correction Functions
% The following two functions pertains to C3 photosynthesis
function [fTv] = temperature_functionC3(Tref,R,T,deltaHa)
% Temperature function
tempfunc1 = (1 - Tref./T);
fTv = exp(deltaHa/(Tref*R).*tempfunc1);
end

function [fHTv] = high_temp_inhibtionC3(Tref,R,T,deltaS,deltaHd)
% High Temperature Inhibition Function
hightempfunc_num = (1+exp((Tref*deltaS-deltaHd)/(Tref*R)));
hightempfunc_deno = (1+exp((deltaS.*T - deltaHd)./(R.*T)));
fHTv = hightempfunc_num ./ hightempfunc_deno;
end

// src/fluxes/biochemical_MD12.m

function biochem_out = biochemical_MD12(leafbio,meteo,~,constants,fV,Q)

%[A,Ci,eta] = biochemical_VCM(Cs,Q,T,eb,O,p,Vcmo,m,Type,Rdparam,stress,Tyear,beta,qLs,NPQs)
% Date:     21 Sep 2012
% Update:   28 Jun 2013 Adaptation for use of Farquhar model of C3 photosynthesis (Farquhar et al 1980)
%           18 Jul 2013 Inclusion of von Caemmerer model of C4 photosynthesis (von Caemmerer 2000, 2013)
%           15 Aug 2013 Modified computation of CO2-limited electron transport in C4 species for consistency with light-limited value
%           22 Oct 2013 Included effect of qLs on Jmax and electron transport; value of kNPQs re-scaled in input as NPQs
%           08 Jan 2019 (CvdT): minor modification to adjust to SCOPE_lite
% Authors:  Federico Magnani, with contributions from Christiaan van der Tol
%
% This function calculates:
%    - CO2 concentration in intercellular spaces (umol/mol == ppmv)
%    - leaf net photosynthesis (umol/m2/s) of C3 or C4 species
%    - fluorescence yield of a leaf (fraction of reference fluorescence yield in dark-adapted and un-stressed leaf)
%
% Usage:
% function [A,Cs,eb,f,rcw] = biochemical(C,Cs,Q,T,ea,eb,O,p,Vcmo,gcparam,Type,tempcor,ra,Tparams,Rdparam,stressfactor,Tyear,beta,qLs,NPQs)
% the function was tested for Matlab 7.2.0.232 (R2006a)
%
% Input (units are important; when not otherwise specified, mol refers to mol C):
% Cs        % [umol/mol]            CO2 concentration at leaf surface
% Q         % [uE/m2/s]             photochemically active radiation absorbed by the leaf
% T         % [oC or K]             leaf temperature
% eb        % [hPa]                 vapour pressure in leaf boundary layer
% O         % [mmol/mol]            ambient O2 concentration
% p         % [Pa]                  air pressure
% Vcmo      % [umol/m2/s]           maximum carboxylation capacity
% m         % [mol/mol]             Ball-Berry coefficient 'm' for stomatal regulation
% Type      % []                    text parameter, either 'C3' for C3 or any other text for C4
% Rdparam   % [mol/mol]             respiration at reference temperature as fraction of Vcmax
% stress    % []                    optional input: stress factor to reduce Vcmax (for example soil moisture, leaf age). Default value = 1 (no stress).
% Tyear     % [oC]                  mean annual temperature
% beta      % []                    fraction of photons partitioned to PSII (0.507 for C3, 0.4 for C4; Yin et al. 2006; Yin and Struik 2012)
% qLs       % []                    fraction of functional reaction centres (Porcar-Castell 2011)
% NPQs      % [s-1]                 rate constant of sustained thermal dissipation, normalized to (kf+kD) (=kNPQs'; Porcar-Castell 2011)
%
% Note: always use the prescribed units. Temperature can be either oC or K
% Note: input can be single numbers, vectors, or n-dimensional matrices
% Note: For consistency reasons, in C4 photosynthesis electron transport rates under CO2-limited conditions are computed by inverting the equation 
%  applied for light-limited conditions(Ubierna et al 2013). A discontinuity would result when computing J from ATP requirements of Vp and Vco, as a 
%  fixed electron transport partitioning is assumed for light-limited conditions

%
% Output:
% A         % [umol/m2/s]           net assimilation rate of the leaves
% Ci        % [umol/mol]            CO2 concentration in intercellular spaces (assumed to be the same as at carboxylation sites in C3 species)
% eta       % []                    amplification factor to be applied to PSII fluorescence yield spectrum
%                                   relative to the dark-adapted, un-stressed yield calculated with either Fluspect or FluorMODleaf

%---------------------------------------------------------------------------------------------------------
%% Start

p       = meteo.p.*1e2;
BallBerrySlope       = leafbio.BallBerrySlope;
BallBerry0          = leafbio.BallBerry0;
O       = meteo.Oa;
Type    = leafbio.Type;
Tyear   = leafbio.Tyear;
beta    = leafbio.beta;
qLs     = leafbio.qLs;
NPQs    = leafbio.kNPQs;
%stress=leafbio.stressfactor;
Cs      = meteo.Cs;
if nargin<6
    Q   = meteo.Q;
end
T       = meteo.T;
eb      = meteo.eb;
Vcmax25    = fV.*leafbio.Vcmax25;
Rdparam = leafbio.RdPerVcmax25;
Q(Q==0) = 1E-9;

%% Global and site-specific constants
R             =  constants.R;                         % [J/K/mol]     universal gas constant

%---------------------------------------------------------------------------------------------------------
%% Unit conversion and computation of environmental variables
T       = T+273.15*(T<100);                           % [K]           convert temperatures to K if not already
RH      = min(1,eb./satvap(T-273.15));                       % []            relative humidity (decimal)
Cs      = Cs .* p .*1E-11;                            % [bar]         1E-6 to convert from ppm to fraction, 1E-5 to convert from Pa to bar
O       = O  .* p .*1E-08;                            % [bar]         1E-3 to convert from mmol/mol to fraction, 1E-5 to convert from Pa to bar

%---------------------------------------------------------------------------------------------------------
%% Define photosynthetic parameters (at reference temperature)
SCOOP     = 2862.;                                    % [mol/mol]     Relative Rubisco specificity for CO2 vs O2 at ref temp (Cousins et al. 2010)
Rdopt     = Rdparam * Vcmax25;                           % [umol/m2/s]   dark respiration at ref temperature from correlation with Vcmax25
switch Type
    case 'C3'                                           % C3 species
        Jmo   =  Vcmax25 * 2.68;                            % [umol/m2/s]   potential e-transport at ref temp from correlation with Vcmax25 (Leuning 1997)
    otherwise                                           % C4 species
        Jmo   =  Vcmax25 * 40/6;                           % [umole-/m2/s] maximum electron transport rate (ratio as in von Caemmerer 2000)
        Vpmo  =  Vcmax25 * 2.33;                             % [umol/m2/s]   maximum PEP carboxylase activity (Yin et al. 2011)
        Vpr   =  80;                                     % [umol/m2/s]   PEP regeneration rate, constant (von Caemmerer 2000)
        gbs   =  (0.0207*Vcmax25+0.4806)*1000.;           % [umol/m2/s]   bundle sheath conductance to CO2 (Yin et al. 2011)
        x     =  0.4;                                     % []            partitioning of electron transport to mesophyll (von Caemmerer 2013)
        alpha =  0;                                      % []            bundle sheath PSII activity (=0 in maize/sorghum; >=0.5 in other cases; von Caemmerer 2000)
end

%---------------------------------------------------------------------------------------------------------
%% Parameters for temperature corrections
TREF         = 25+273.15;                             % [K]            reference temperature for photosynthetic processes

HARD         = 46.39;                                 % [kJ/mol]       activation energy of Rd
CRD          = 1000.*HARD/(R*TREF);                   % []             scaling factor in RD response to temperature

HAGSTAR      = 37.83;                                 % [kJ/mol]       activation energy of Gamma_star
CGSTAR       = 1000.*HAGSTAR/(R*TREF);                % []             scaling factor in GSTAR response to temperature

switch Type
    case 'C3'                                           % C3 species
        HAJ     = 49.88;                                 % [kJ/mol]       activation energy of Jm (Kattge & Knorr 2007)
        HDJ     = 200;                                   % [kJ/mol]       deactivation energy of Jm (Kattge & Knorr 2007)
        DELTASJ = (-0.75*Tyear+660)/1000;                % [kJ/mol/K]     entropy term for J  (Kattge and Knorr 2007)
        
        HAVCM   = 71.51;                                 % [kJ/mol]       activation energy of Vcm (Kattge and Knorr 2007)
        HDVC    = 200;                                   % [kJ/mol]       deactivation energy of Vcm (Kattge & Knorr 2007)
        DELTASVC= (-1.07*Tyear+668)/1000;                % [kJ/mol/K]     entropy term for Vcmax (Kattge and Knorr 2007)
        
        KCOP    = 404.9;                                 % [umol/mol]     Michaelis-Menten constant for CO2 at ref temp (Bernacchi et al 2001)
        HAKC    = 79.43;                                 % [kJ/mol]       activation energy of Kc (Bernacchi et al 2001)
        
        KOOP    = 278.4;                                 % [mmol/mol]     Michaelis-Menten constant for O2  at ref temp (Bernacchi et al 2001)
        HAKO    = 36.38;                                 % [kJ/mol]       activation energy of Ko (Bernacchi et al 2001)
        
    otherwise                                           % C4 species (values can be different as noted by von Caemmerer 2000)
        HAJ    = 77.9;                                   % [kJ/mol]       activation energy of Jm  (Massad et al 2007)
        HDJ     = 191.9;                                 % [kJ/mol]       deactivation energy of Jm (Massad et al 2007)
        DELTASJ = 0.627;                                 % [kJ/mol/K]     entropy term for Jm (Massad et al 2007). No data available on acclimation to temperature.
        
        HAVCM   = 67.29;                                 % [kJ/mol]       activation energy of Vcm (Massad et al 2007)
        HDVC    = 144.57;                                % [kJ/mol]       deactivation energy of Vcm (Massad et al 2007)
        DELTASVC= 0.472;                                 % [kJ/mol/K]     entropy term for Vcm (Massad et al 2007). No data available on acclimation to temperature.
        
        HAVPM   = 70.37;                                 % [kJ/mol]       activation energy of Vpm  (Massad et al 2007)
        HDVP    = 117.93;                                % [kJ/mol]       deactivation energy of Vpm (Massad et al 2007)
        DELTASVP= 0.376;                                 % [kJ/mol/K]     entropy term for Vpm (Massad et al 2007). No data available on acclimation to temperature.
        
        KCOP    = 944.;                                  % [umol/mol]     Michaelis-Menten constant for CO2 at ref temp (Chen et al 1994; Massad et al 2007)
        Q10KC   = 2.1;                                   % []             Q10 for temperature response of Kc (Chen et al 1994; Massad et al 2007)
        
        KOOP    = 633.;                                  % [mmol/mol]     Michaelis-Menten constant for O2 at ref temp (Chen et al 1994; Massad et al 2007)
        Q10KO   = 1.2;                                   % []             Q10 for temperature response of Ko (Chen et al 1994; Massad et al 2007)
        
        KPOP    = 82.;                                   % [umol/mol]     Michaelis-Menten constant of PEP carboxylase at ref temp (Chen et al 1994; Massad et al 2007)
        Q10KP   = 2.1;                                   % []             Q10 for temperature response of Kp (Chen et al 1994; Massad et al 2007)
        
end


%---------------------------------------------------------------------------------------------------------
%% Corrections for effects of temperature and non-stomatal limitations
dum1   = R./1000.*T;                                  % [kJ/mol]
dum2   = R./1000.*TREF;                               % [kJ/mol]

Rd     = Rdopt.*exp(CRD-HARD./dum1);                  % [umol/m2/s]    mitochondrial respiration rates adjusted for temperature (Bernacchi et al. 2001)
SCO    = SCOOP./exp(CGSTAR-HAGSTAR./dum1);            % []             Rubisco specificity for CO2 adjusted for temperature (Bernacchi et al. 2001)

Jmax   = Jmo .* exp(HAJ.*(T-TREF)./(TREF*dum1));
Jmax   = Jmax.*(1.+exp((TREF*DELTASJ-HDJ)./dum2));
Jmax   = Jmax./(1.+exp((T.*DELTASJ-HDJ)./dum1));     % [umol e-/m2/s] max electron transport rate at leaf temperature (Kattge and Knorr 2007; Massad et al. 2007)

Vcmax  = Vcmax25 .* exp(HAVCM.*(T-TREF)./(TREF*dum1));
Vcmax  = Vcmax.*(1+exp((TREF*DELTASVC-HDVC)/dum2));
Vcmax  = Vcmax./(1+exp((T.*DELTASVC-HDVC)./dum1));    % [umol/m2/s]    max carboxylation rate at leaf temperature (Kattge and Knorr 2007; Massad et al. 2007)

switch Type
    case 'C3'                                           % C3 species
        CKC    = 1000.*HAKC/(R*TREF);                     % []             scaling factor in KC response to temperature
        Kc     = KCOP.*exp(CKC-HAKC./dum1).*1e-11.*p;     % [bar]          Michaelis constant of carboxylation adjusted for temperature (Bernacchi et al. 2001)
        
        CKO    = 1000.*HAKO/(R*TREF);                     % []             scaling factor in KO response to temperature
        Ko     = KOOP.*exp(CKO-HAKO./dum1).*1e-8.*p;      % [bar]          Michaelis constant of oxygenation adjusted for temperature (Bernacchi et al. 2001)
        
    otherwise                                           % C4 species
        Vpmax  = Vpmo .* exp(HAVPM.*(T-TREF)./(TREF*dum1));
        Vpmax  = Vpmax.*(1+exp((TREF*DELTASVP-HDVP)/dum2));
        Vpmax  = Vpmax./(1+exp((T.*DELTASVP-HDVP)./dum1));% [umol/m2/s]    max carboxylation rate at leaf temperature (Massad et al. 2007)
        
        Kc     = KCOP.*Q10KC .^ ((T-TREF)/10.)*1e-11*p;    % [bar]          Michaelis constant of carboxylation temperature corrected (Chen et al 1994; Massad et al 2007)
        
        Ko     = KOOP.*Q10KO .^ ((T-TREF)/10.)*1e-8*p;     % [bar]          Michaelis constant of oxygenation  temperature corrected (Chen et al 1994; Massad et al 2007)
        
        Kp     = KPOP.*Q10KP .^ ((T-TREF)/10.)*1e-11*p;    % [bar]          Michaelis constant of PEP carboxyl temperature corrected (Chen et al 1994; Massad et al 2007)
        
end

%---------------------------------------------------------------------------------------------------------
%% Define electron transport and fluorescence parameters
kf        = 3.E7;                                    % [s-1]         rate constant for fluorescence
kD        = 1.E8;                                   % [s-1]         rate constant for thermal deactivation at Fm
kd        = 1.95E8;                                    % [s-1]         rate constant of energy dissipation in closed RCs (for theta=0.7 under un-stressed conditions)  
po0max    = 0.88;                                     % [mol e-/E]    maximum PSII quantum yield, dark-acclimated in the absence of stress (Pfundel 1998)
kPSII     = (kD+kf) * po0max/(1.-po0max);             % [s-1]         rate constant for photochemisty (Genty et al. 1989)
fo0       = kf./(kf+kPSII+kD);                        % [E/E]         reference dark-adapted PSII fluorescence yield under un-stressed conditions

kps       = kPSII * qLs;                              % [s-1]         rate constant for photochemisty under stressed conditions (Porcar-Castell 2011)
kNPQs     = NPQs * (kf+kD);                           % [s-1]         rate constant of sustained thermal dissipation (Porcar-Castell 2011)
kds       = kd * qLs;
kDs       = kD + kNPQs;
Jms       = Jmax * qLs;                               % [umol e-/m2/s] potential e-transport rate reduced for PSII photodamage
po0       = kps ./(kps+kf+kDs);                       % [mol e-/E]    maximum PSII quantum yield, dark-acclimated in the presence of stress
THETA     = (kps-kds)./(kps+kf+kDs);                  % []            convexity factor in J response to PAR

%---------------------------------------------------------------------------------------------------------
%% Calculation of electron transport rate
Q2     = beta * Q * po0;
J      = (Q2+Jms-sqrt((Q2+Jms).^2-4*THETA.*Q2.*Jms))./(2*THETA); % [umol e-/m2/s]    electron transport rate under light-limiting conditions

%---------------------------------------------------------------------------------------------------------
%% Calculation of net photosynthesis
switch Type
    case 'C3'
        minCi = 0.3;
    otherwise
        minCi = 0.1;
end

Ci = BallBerry(Cs, RH, [], BallBerrySlope, BallBerry0, minCi);

switch Type
    case 'C3'                                           % C3 species, based on Farquhar model (Farquhar et al. 1980)
        GSTAR = 0.5*O./SCO;                             % [bar]             CO2 compensation point in the absence of mitochondrial respiration
        
        Cc  = Ci;                                        % [bar]             CO2 concentration at carboxylation sites (neglecting mesophyll resistance)
        
        Wc  = Vcmax .* Cc ./ (Cc + Kc .* (1+O./Ko));     % [umol/m2/s]       RuBP-limited carboxylation
        Wj  = J.*Cc ./ (4.5*Cc + 10.5*GSTAR);            % [umol/m2/s]       electr transp-limited carboxyl
        
        W   = min(Wc,Wj);                                % [umol/m2/s]       carboxylation rate
        Ag  = (1 - GSTAR./Cc) .*W;                       % [umol/m2/s]       gross photosynthesis rate
        A   = Ag - Rd;                                   % [umol/m2/s]       net photosynthesis rate
        Ja  = J.*W ./Wj;                                 % [umole-/m2/s]     actual linear electron transport rate
        
    otherwise                                           % C4 species, based on von Caemmerer model (von Caemmerer 2000)
        %Ci    =  max(9.9e-6*(p*1e-5),Cs.*(1-1.6./(m.*RH*stress)));
        % [bar]             intercellular CO2 concentration from Ball-Berry model (Ball et al. 1987)
        minCi = 0.1;
        Ci = BallBerry(Cs, RH, [], BallBerrySlope, 0, minCi);
        
        Cm    =  Ci;                                     % [bar]             mesophyll CO2 concentration (neglecting mesophyll resistance)
        Rs    =  0.5 .* Rd;                               % [umol/m2/s]       bundle sheath mitochondrial respiration (von Caemmerer 2000)
        Rm    =  Rs;                                     % [umol/m2/s]       mesophyll mitochondrial respiration
        gam   =  0.5./SCO;                               % []                half the reciprocal of Rubisco specificity for CO2
        
        Vpc   = Vpmax .* Cm./(Cm+Kp);                     % [umol/m2/s]       PEP carboxylation rate under limiting CO2 (saturating PEP)
        Vp    = min(Vpc,Vpr);                            % [umol/m2/s]       PEP carboxylation rate
        
        % Complete model proposed by von Caemmerer (2000)
        dum1  =  alpha/0.047;                           % dummy variables, to reduce computation time
        dum2  =  Kc./Ko;
        dum3  =  Vp-Rm+gbs.*Cm;
        dum4  =  Vcmax-Rd;
        dum5  =  gbs.*Kc.*(1+O./Ko);
        dum6  =  gam.*Vcmax;
        dum7  =  x*J./2. - Rm + gbs.*Cm;
        dum8  =  (1.-x).*J./3.;
        dum9  =  dum8 - Rd;
        dum10 =  dum8 + Rd * 7/3;
        
        a     =  1. - dum1 .* dum2;
        b     =  -(dum3+dum4+dum5+dum1.*(dum6+Rd.*dum2));
        c     =  dum4.*dum3-dum6.*gbs*O+Rd.*dum5;
        Ac    =  (-b - sqrt(b.^2-4.*a.*c))./(2.*a);           % [umol/m2/s]       CO2-limited net photosynthesis
        
        a     =  1.- 7./3.*gam.*dum1;
        b     =  -(dum7+dum9 + gbs.*gam.*O.*7./3. + dum1.*gam.*dum10);
        c     =  dum7.*dum9 - gbs.*gam.*O.*dum10;
        Aj    =  (-b - sqrt(b.^2-4.*a.*c))./(2.*a);           % [umol/m2/s]       light-limited net photosynthesis (assuming that an obligatory Q cycle operates)
        
        A     =  min(Ac,Aj);                             % [umol/m2/s]       net photosynthesis
        
        Ja    =  J;                                      % [umole-/m2/s]     actual electron transport rate, CO2-limited
               
        if any(A==Ac) %IPL 03/09/2013
            
            ind=A==Ac;
            a(ind)   =  x.*(1-x)./6./A(ind);
            b(ind)   =  (1-x)/3.*(gbs(ind)./A(ind).*(Cm(ind)-Rm(ind)./gbs(ind)-gam(ind).*O)-1-alpha.*gam(ind)./0.047)-x./2.*(1.+Rd(ind)./A(ind));
            c(ind)   =  (1+Rd(ind)./A(ind)).*(Rm(ind)-gbs(ind).*Cm(ind)-7.*gbs(ind).*gam(ind).*O./3)+(Rd(ind)+A(ind)).*(1-7.*alpha.*gam(ind)./3./0.047);
            Ja(ind)  =  (-b(ind) + sqrt(b(ind).^2-4.*a(ind).*c(ind)))./(2.*a(ind));            % [umole-/m2/s]     actual electron transport rate, CO2-limited
            
        end
end

%---------------------------------------------------------------------------------------------------------
%% Calculation of PSII quantum yield and fluorescence
ps     = Ja ./(beta.*Q);                            % [mol e-/E]    PSII photochemical quantum yield
[fs]   = MD12(ps,Ja,Jms,kps,kf,kds,kDs);            % [E/E]         PSII fluorescence yield
eta    = fs./fo0;                                   % []            scaled PSII fluorescence yield

%% JP add
rhoa        = 1.2047;                               % [kg m-3]       specific mass of air
Mair        = 28.96;                                % [g mol-1]      molecular mass of dry air

%rcw         = 0.625*(Cs-Ci)./A *rhoa/Mair*1E3    * 1e6 ./ p .* 1E5;
ppm2bar =  1e-6 .* (p .*1E-5);
gs = 1.6 * A* ppm2bar./ (Cs-Ci);
rcw      =  (rhoa./(Mair*1E-3))./gs;
rcw(A<=0 & rcw~=0)   = 0.625*1E6;

%% convert back to ppm
Ci          = Ci*1e6 ./ p .* 1E5;

%%
biochem_out.A = A;
biochem_out.Ci = Ci;
biochem_out.ps = ps;
biochem_out.eta = eta;
biochem_out.fs  = fs;
biochem_out.rcw = rcw;
biochem_out.qE  = rcw*NaN; % dummy output, to be consistent with SCOPE
biochem_out.Kn  = NPQs + 0*rcw; %
biochem_out.Phi_N  = kNPQs./(kNPQs +kD+kf+kps)+ 0*rcw;
biochem_out.Ja      = Ja;

return;

end
%%% end of function biochemical


%---------------------------------------------------------------------------------------------------------
%% MD12 algorithm for the computation of fluorescence yield

function [fs] = MD12(ps,Ja,Jms,kps,kf,kds,kDs)

fs1    = ps .* (kf./kps) ./ (1. - Ja./Jms);         % [E/E]   PSII fluorescence yield under CO2-limited conditions

par1   = kps./(kps-kds);                            % [E/E]   empirical parameter in the relationship under light-limited conditions
par2   = par1.* (kf+kDs+kds)./kf;                  % [E/E]   empirical parameter in the relationship under light-limited conditions
fs2    = (par1-ps)./par2;                           % [E/E]   PSII fluorescence yield under light-limited conditions

fs     = min(fs1,fs2);                              % [E/E]   PSII fluorescence yield
end



%% Ball Berry Model
function [Ci, gs] = BallBerry(Cs, RH, A, BallBerrySlope, BallBerry0, minCi, Ci_input)
%  Cs  : CO2 at leaf surface
%  RH  : relative humidity
%  A   : Net assimilation in 'same units of CO2 as Cs'/m2/s
% BallBerrySlope, BallBerry0, 
% minCi : minimum Ci as a fraction of Cs (in case RH is very low?)
% Ci_input : will calculate gs if A is specified.
if nargin > 6 && ~isempty(Ci_input)
    % Ci is given: try and compute gs
    Ci = Ci_input;
    gs = [];
    if ~isempty(A) && nargout > 1
        gs = gsFun(Cs, RH, A, BallBerrySlope, BallBerry0);
    end
elseif all(BallBerry0 == 0) || isempty(A)
    % EXPLANATION:   *at equilibrium* CO2_in = CO2_out => A = gs(Cs - Ci) [1]
    %  so Ci = Cs - A/gs (at equilibrium)                                 [2]
    %  Ball-Berry suggest: gs = m (A RH)/Cs + b   (also at equilib., see Leuning 1990)
    %  if b = 0 we can rearrange B-B for the second term in [2]:  A/gs = Cs/(m RH)
    %  Substituting into [2]
    %  Ci = Cs - Cs/(m RH) = Cs ( 1- 1/(m RH)  [ the 1.6 converts from CO2- to H2O-diffusion ]
    Ci      = max(minCi .* Cs,  Cs.*(1-1.6./(BallBerrySlope .* RH)));
    gs      = [];
else
    %  if b > 0  Ci = Cs( 1 - 1/(m RH + b Cs/A) )
    % if we use Leuning 1990, Ci = Cs - (Cs - Gamma)/(m RH + b(Cs - Gamma)/A)  [see def of Gamma, above]
    % note: the original B-B units are A: umol/m2/s, ci ppm (umol/mol), RH (unitless)
    %   Cs input was ppm but was multiplied by ppm2bar above, so multiply A by ppm2bar to put them on the same scale.
    %  don't let gs go below its minimum value (i.e. when A goes negative)
    gs = gsFun(Cs, RH, A, BallBerrySlope, BallBerry0);
    Ci = max(minCi .* Cs,  Cs - 1.6 * A./gs) ;
end

end % function

function gs = gsFun(Cs, RH, A, BallBerrySlope, BallBerry0)
% add in a bit just to avoid div zero. 1 ppm = 1e-6 (note since A < 0 if Cs ==0, it gives a small gs rather than maximal gs
gs = max(BallBerry0,  BallBerrySlope.* A .* RH ./ (Cs+1e-9)  + BallBerry0);
% clean it up:
%gs( Cs == 0 ) = would need to be max gs here;  % eliminate infinities
gs( isnan(Cs) ) = NaN;  % max(NaN, X) = X  (MATLAB 2013b) so fix it here
end




% Sources:
%  Ball J. T., I. E. Woodrow and J. A. Berry. (1987) A model predicting stomatal conductance and its contribution to the control of photosynthesis
%    under different environmental conditions. In: Progress in Photosynthesis Research (Ed. J. Biggens), p. 221-224, The Netherlands:Martinus Nijhoff.
%  Bernacchi C.J., E.L. Singsaas, C. Pimentel, A.R. Portis and S.P. Long (2001) Improved temperature response functions for models of Rubisco-limited
%    photosynthesis. Plant Cell Envir 24:253-259.
%  Bernacchi C.J., C. Pimentel and S.P. Long (2003) In vivo temperature response functions of parameters required to model RuBP-limited photosynthesis.
%    Plant Cell Envir 26 (9):1419-1430.
%  Chen D.X., M.B. Coughenour, A.K. Knapp, and C.E. Owensby (1994) Mathematical simulation of C4 grass photosynthesis in ambient and elevated CO2.
%    Ecol.Model. 73:63-80, 1994.
%  Cousins A.B., O. Ghannoum, S. von Caemmerer, and M.R. Badger (2010) Simultaneous determination of Rubisco carboxylase and oxygenase kinetic parameters
%    in Triticum aestivum and Zea mays using membrane inlet mass spectrometry. Plant Cell Envir 33:444-452.
%  Farquhar G.D., S. von Caemmerer and J.A. Berry (1980) A biochemical model of photosynthetic CO2 assimilation in leaves of C3 species. Planta 149:78-90.
%  Genty B., J.-M. Briantais and N. R. Baker (1989) The relationship between quantum yield of photosynthetic electron transport and quenching of
%    chlorophyll fluorescence. Biochimica et Biophysica Acta 990:87-92.
%  Kattge  J. and W. Knorr (2007) Temperature acclimation in a biochemical model of photosynthesis: a reanalysis of data from 36 species.
%    Plant Cell Envir 30:1176-1190.
%  Leuning R. (1997) Scaling to a common temperature improves the correlation between the photosynthesis parameters Jmax and Vcmax.
%    J.Exp.Bot. 48 (307):345-347.
%  Massad R.S., A. Tuzet and O. Bethenod (2007) The effect of temperature on C4-type leaf photosynthesis parameters. Plant Cell Envir 30:1191-1204.
%  Pfundel E. (1998) Estimating the contribution of Photosystem I to total leaf chlorophyll fluorescence. Photosynthesis Research 56:185-195.
%  Porcar-Castell A. (2011) A high-resolution portrait of the annual dynamics of photochemical and non-photochemical quenching in needles of  Pinus sylvestris.
%    Physiol.Plant. 143:139-153.
%  von Caemmerer S. (2000) Biochemical Models of Leaf Photosynthesis, Canberra:CSIRO Publishing.
%  von Caemmerer S. (2013) Steady-state models of photosynthesis. Plant Cell Envir, in press.
%  Yin X., Z. Sun, P.C. Struik, P.E.L. van der Putten, W. van Ieperen and J. Harbinson (2011) Using a biochemical C4 photosynthesis model and combined
%    gas exchange and chlorophyll fluorescence measurements to estimate bundle-sheath conductance of maize leaves differing in age and nitrogen content.
%    Plant Cell Envir 34:2183-2199.
%

// src/fluxes/ebal.m

function [iter,rad,thermal,soil,bcu,bch,fluxes,resist_out,meteo]             ...
    = ebal(constants,options,rad,gap,  ...
    meteo,soil,canopy,leafbio,k,xyt,integr)
% function ebal.m calculates the energy balance of a vegetated surface
%
% authors:      Christiaan van der Tol (c.vandertol@utwente.nl)
%               Joris Timmermans
% date          26 Nov 2007 (CvdT)
% updates       29 Jan 2008 (JT & CvdT)     converted into a function
%               11 Feb 2008 (JT & CvdT)     improved soil heat flux and temperature calculation
%               14 Feb 2008 (JT)            changed h in to hc (as h=Avogadro`s constant)
%               31 Jul 2008 (CvdT)          Included Pntot in output
%               19 Sep 2008 (CvdT)          Converted F0 and F1 from units per aPAR into units per iPAR
%               07 Nov 2008 (CvdT)          Changed layout
%               18 Sep 2012 (CvdT)          Changed Oc, Cc, ec
%                  Feb 2012 (WV)            introduced structures for variables
%                  Sep 2013 (JV, CvT)       introduced additional biochemical model
%               10 Dec 2019 (CvdT)          made a light version (layer
%                                           averaged fluxes)
% parent: SCOPE.m (script)
% uses:
%       RTMt_sb.m, RTMt_planck.m (optional), RTMf.m (optional)
%       resistances.m
%       heatfluxes.m
%       biochemical.m
%       soil_respiration.m
%
% Table of contents of the function
%
%   1. Initialisations for the iteration loop
%           intial values are attributed to variables
%   2. Energy balance iteration loop
%           iteration between thermal RTM and surface fluxes
%   3. Write warnings whenever the energy balance did not close
%   4. Calculate vertical profiles (optional)
%   5. Calculate spectrally integrated energy, water and CO2 fluxes
%
% The energy balance iteration loop works as follows:
%
% RTMo              More or less the classic SAIL model for Radiative
%                   Transfer of sun and sky light (no emission by the vegetation)
% While continue	Here an iteration loop starts to close the energy
%                   balance, i.e. to match the micro-meteorological model
%                   and the radiative transfer model
% 	RTMt_sb         A numerical Radiative Transfer Model for thermal
%                   radiation emitted by the vegetation
% 	resistances     Calculates aerodynamic and boundary layer resistances
%                   of vegetation and soil (the micro-meteorological model)
% 	biochemical     Calculates photosynthesis, fluorescence and stomatal
%                   resistance of leaves (or biochemical_MD12: alternative)
% 	heatfluxes      Calculates sensible and latent heat flux of soil and
%                   vegetation
%                   Next soil heat flux is calculated, the energy balance
%                   is evaluated, and soil and leaf temperatures adjusted
%                   to force energy balance closure
% end {while continue}
%
% meanleaf          Integrates the fluxes over all leaf inclinations
%                   azimuth angles and layers, integrates over the spectrum
%
% usage:
%[iter,fluxes,rad,profiles,thermal]             ...
%         = ebal(iter,options,spectral,rad,gap,leafopt,  ...
%                angles,meteo,soil,canopy,leafbio)
%
% The input and output are structures. These structures are further
% specified in a readme file.
%
% Input:
%
%   iter        numerical parameters used in the iteration for energy balance closure
%   options     calculation options
%   spectral    spectral resolutions and wavelengths
%   rad         incident radiation
%   gap         probabilities of direct light penetration and viewing
%   leafopt     leaf optical properties
%   angles      viewing and observation angles
%   soil        soil properties
%   canopy      canopy properties
%   leafbio     leaf biochemical parameters
%
% Output:
%
%   iter        numerical parameters used in the iteration for energy balance closure
%   fluxes      energy balance, turbulent, and CO2 fluxes
%   rad         radiation spectra
%   thermal     temperatures, aerodynamic resistances and friction velocity
%   bcu, bch    leaf biochemical outputs for sunlit and shaded leaves,
%               respectively

%% 1. initialisations and other preparations for the iteration loop
% parameters for the closure loop
counter     = 0;                % iteration counter of ebal
maxit       = 100;              % maximum number of iterations
maxEBer     = 1;                % maximum energy balance error (any leaf) [Wm-2]
Wc          = 1;                % update step (1 is nominal, [0,1] possible)
CONT        = 1;                % boolean indicating whether iteration continues

% constants
MH2O        = constants.MH2O;
Mair        = constants.Mair;
rhoa        = constants.rhoa;
cp          = constants.cp;
sigmaSB     = constants.sigmaSB;

% input preparation
nl          = canopy.nlayers;   
GAM         = soil.GAM;
Ps          = gap.Ps;
kV          = canopy.kV;
xl          = canopy.xl;
LAI         = canopy.LAI;
rss         = soil.rss;

% functions for saturated vapour pressure 
es_fun      = @(T)6.107*10.^(7.5.*T./(237.3+T));
s_fun       = @(es, T) es*2.3026*7.5*237.3./(237.3+T).^2;

SoilHeatMethod = options.soil_heat_method;
if ~(options.simulation==1), SoilHeatMethod = 2; end
if SoilHeatMethod < 2
    if k > 1
        Deltat          = (datenum(xyt.t(k))-datenum(xyt.t(k-1)))*86400;           %           Duration of the time interval (s)
    else
        Deltat          = 1/48*86400;
    end   
    x 		= [1:12;1:12]'*Deltat;
    Tsold   = soil.Tsold;
end

% meteo
Ta          = meteo.Ta;
ea          = meteo.ea;
Ca          = meteo.Ca;
p           = meteo.p;
Rnuc        = rad.Rnuc;
ech         = ea*ones(nl,1);          % Leaf boundary vapour pressure (shaded/sunlit leaves)
Cch         = Ca*ones(nl,1);
ecu         = ea+0*Rnuc;
Ccu         = Ca+0*Rnuc;          % Leaf boundary CO2 (shaded/sunlit leaves)

% other preparations
e_to_q      = MH2O/Mair./p;             % Conversion of vapour pressure [Pa] to absolute humidity [kg kg-1]
Fc          = Ps(1:end-1);
Fs          = [1-Ps(end),Ps(end)];      % Matrix containing values for 1-Ps and Ps of soil
%Fc          = (1-Ps(1:end-1))'/nl;      % Matrix containing values for Ps of canopy
fV          = exp(kV*xl(1:end-1));      % Vertical profile of Vcmax

% initial values for the loop
Ts          = (Ta+3)*ones(2,1);         % soil temperature (+3 for a head start of the iteration) 
Tch         = (Ta+.1)*ones(nl,1);       % leaf temperature (shaded leaves)
Tcu         = (Ta+.3)*ones(size(Rnuc)); % leaf tempeFrature (sunlit leaves)
meteo.L     = -1E6;                     % Monin-Obukhov length
[meteo_h,meteo_u]  = deal(meteo);

% this for is the exponential decline of Vcmax25. If 'lite' the dimensions
% are [nl], otherwise [13,36,nl]: 
if size(Rnuc,2)>1
    fVu      = ones(13,36,nl);
    for i = 1:nl
        fVu(:,:,i) = fV(i);
    end
else
    fVu = fV;
end

%% 2.1 Energy balance iteration loop
%Energy balance loop (Energy balance and radiative transfer)

while CONT                          % while energy balance does not close
    % 2.1. Net radiation of the components
    % Thermal radiative transfer model for vegetation emission (with Stefan-Boltzman's equation)
    rad     = RTMt_sb(constants,rad,soil,leafbio,canopy,gap,Tcu,Tch,Ts(2),Ts(1),0);
    Rnhc    = rad.Rnhc + rad.Rnhct;     % Canopy (shaded) net radiation
    Rnuc    = rad.Rnuc + rad.Rnuct;     % Canopy (sunlit) net radiation
    Rnhs    = rad.Rnhs+rad.Rnhst;       % Soil   (sun+sh) net radiation
    Rnus    = rad.Rnus+rad.Rnust;
    Rns     = [Rnhs Rnus]';
    
    % 2.3. Biochemical processes
    meteo_h.T       = Tch;
    meteo_h.eb      = ech;
    meteo_h.Cs      = Cch;
    meteo_h.Q       = rad.Pnh_Cab;
    meteo_u.T       = Tcu;
    meteo_u.eb      = ecu;
    meteo_u.Cs      = Ccu;
    meteo_u.Q       = rad.Pnu_Cab;
    if options.Fluorescence_model == 1
        b       = @biochemical_MD12;
    else
        b       = @biochemical;
    end
    bch     = b(leafbio,meteo_h,options,constants,fV);
    bcu     = b(leafbio,meteo_u,options,constants,fVu);
     
    % Aerodynamic roughness
    % calculate friction velocity [m s-1] and aerodynamic resistances [s m-1]  
    [resist_out]  = resistances(constants,soil,canopy,meteo);
    meteo.ustar = resist_out.ustar;
    raa     = resist_out.raa;
    rawc    = resist_out.rawc;
    raws    = resist_out.raws;  
    rac     = (LAI+1)*(raa+rawc);
    ras     = (LAI+1)*(raa+raws);
    
    % Fluxes (latent heat flux (lE), sensible heat flux (H) and soil heat flux G
    % in analogy to Ohm's law, for canopy (c) and soil (s). All in units of [W m-2]
    [lEch,Hch,ech,Cch,lambdah,sh]     = heatfluxes(rac,bch.rcw,Tch,ea,Ta,e_to_q,Ca,bch.Ci,constants, es_fun, s_fun);
    [lEcu,Hcu,ecu,Ccu,lambdau,su]     = heatfluxes(rac,bcu.rcw,Tcu,ea,Ta,e_to_q,Ca,bcu.Ci,constants, es_fun, s_fun);
    [lEs,Hs,~,~,lambdas,ss]           = heatfluxes(ras,rss ,Ts ,ea,Ta,e_to_q,Ca,Ca,constants, es_fun, s_fun);
    
    % integration over the layers and sunlit and shaded fractions
    Hstot   = Fs*Hs;
    %Hctot   = LAI*(meanleaf(canopy,Hcu,integr,Ps(1:end-1)) + meanleaf(canopy,Hch,'layers',1-Ps(1:end-1))  ); 
    Hctot   = aggregator(LAI,Hcu, Hch, Ps(1:end-1), canopy,integr);
    
    Htot    = Hstot + Hctot;
    if options.MoninObukhov
        meteo.L     = Monin_Obukhov(constants,meteo,Htot);     
    end
    
    % ground heat flux
    if SoilHeatMethod == 2
        G = 0.35*Rns;
        dG = 4*(1-soil.rs_thermal)*sigmaSB*(Ts+273.15).^3 * 0.35;
    elseif SoilHeatMethod == 3  % daily average flux
        G = 0*Rns;
        dG = 0*Rns;
    else
        G = GAM/sqrt(pi) * 2* sum(([Ts'; Tsold(1:end-1,:)] - Tsold)/Deltat .* (sqrt(x) - sqrt(x-Deltat)));
        G = G';
        dG = GAM/sqrt(pi) * 2* ((sqrt(x(1)) - sqrt(x(1)-Deltat)))/Deltat * ones(2,1);
    end
    
    % energy balance errors, continue criterion and iteration counter
    EBerch  = Rnhc -lEch -Hch;
    EBercu  = Rnuc -lEcu -Hcu;
    EBers   = Rns  -lEs  -Hs - G;
    
    counter     = counter+1;                   %        Number of iterations
    maxEBercu   = max(max(max(abs(EBercu))));
    maxEBerch   = max(abs(EBerch));
    maxEBers    = max(abs(EBers));
    
    CONT        = ( maxEBercu >   maxEBer    |...
        maxEBerch >   maxEBer     |...
        maxEBers  >   maxEBer)    &...
        counter   <   maxit+1;%        Continue iteration?
    if ~CONT
        if any(isnan([maxEBercu, maxEBerch, maxEBers]))
            fprintf('WARNING: NaN in fluxes, counter = %i\n', counter)
        end
        break
    end
    if counter==10, Wc = 0.8;  end
    if counter==20; Wc = 0.6;  end

    % if counter>99, plot(EBercu(:)), hold on, end
    % 2.7. New estimates of soil (s) and leaf (c) temperatures, shaded (h) and sunlit (1)
    Tch         = Tch + Wc*EBerch./((rhoa*cp)./rac + rhoa*lambdah*e_to_q.*sh./(rac+bch.rcw)+ 4*leafbio.emis*sigmaSB*(Tch+273.15).^3);
    Tcu         = Tcu + Wc*EBercu./((rhoa*cp)./rac + rhoa*lambdau*e_to_q.*su./(rac+bcu.rcw)+ 4*leafbio.emis*sigmaSB*(Tcu+273.15).^3);
    Ts          = Ts + Wc*EBers./(rhoa*cp./ras + rhoa*lambdas*e_to_q.*ss/(ras+rss)+ 4*(1-soil.rs_thermal)*sigmaSB*(Ts+273.15).^3 + dG);
    Tch(abs(Tch)>100) = Ta;
    Tcu(abs(Tcu)>100) = Ta;
end

%% 2.2 emmissivity calculation
rad     = RTMt_sb(constants,rad,soil,leafbio,canopy,gap,Tcu,Tch,Ts(2),Ts(1),0);
[blackleaf.tau_thermal,blackleaf.rho_thermal,blacksoil.rs_thermal] = deal(0);
rad0    = RTMt_sb(constants,rad,blacksoil,blackleaf,canopy,gap,Tcu,Tch,Ts(2),Ts(1),0);
rad.canopyemis = rad.Eoutte./rad0.Eoutte;

%% 3. Print warnings whenever the energy balance could not be solved
if counter>=maxit
    fprintf('WARNING: maximum number of iteratations exceeded\n');
    fprintf('Maximum energy balance error sunlit vegetation = %4.2f W m-2\n' ,maxEBercu);
    fprintf('Maximum energy balance error shaded vegetation = %4.2f W m-2\n' ,maxEBerch);
    fprintf('Energy balance error soil              = %4.2f W m-2\n' ,maxEBers);
    fprintf('Mean error sunlit vegetation = %4.2f W m-2\n' , mean(EBercu(:)));
end

%% 4. some more outputs
iter.counter    = counter;

thermal.Tcu     = Tcu;
thermal.Tch     = Tch;
thermal.Tsu     = Ts(2);
thermal.Tsh     = Ts(1);

fluxes.Rnctot = aggregator(LAI,Rnuc, Rnhc, Fc, canopy,integr);     % net radiation leaves
fluxes.lEctot = aggregator(LAI,lEcu, lEch, Fc, canopy,integr);     % latent heat leaves
fluxes.Hctot  = aggregator(LAI,Hcu, Hch, Fc,  canopy,integr);       % sensible heat leaves
fluxes.Actot  = aggregator(LAI,bcu.A, bch.A, Fc, canopy,integr);   % photosynthesis leaves
fluxes.Tcave  = aggregator(1,Tcu, Tch, Fc,  canopy,integr);         % mean leaf temperature
fluxes.Rnstot = Fs*Rns;           % Net radiation soil
fluxes.lEstot = Fs*lEs;           % Latent heat soil
fluxes.Hstot  = Fs*Hs;            % Sensible heat soil
fluxes.Gtot   = Fs*G;             % Soil heat flux
fluxes.Tsave  = Fs*Ts;            % Soil temperature
% fluxes.Resp   = Fs*equations.soil_respiration(Ts); %  Soil respiration = 0
fluxes.Rntot = fluxes.Rnctot + fluxes.Rnstot;
fluxes.lEtot = fluxes.lEctot + fluxes.lEstot;
fluxes.Htot = fluxes.Hctot + fluxes.Hstot;

resist_out.rss = rss; % this is simply a copy of the input rss

%% update soil temperatures history
if SoilHeatMethod < 2
    Tsold(2:end,:) = soil.Tsold(1:end-1,:);
    Tsold(1,:) 	= Ts(:);
    if isnan(Ts)
        Tsold(1,:) = Tsold(2,:);
    end
    soil.Tsold = Tsold;
end
return

function flux_tot = aggregator(LAI,sunlit_flux, shaded_flux, Fs, canopy,integr)
flux_tot = LAI*(meanleaf(canopy,sunlit_flux,integr,Fs) + meanleaf(canopy,shaded_flux,'layers',1-Fs));
return

// src/fluxes/ebal_bigleaf.m

function [iter,rad,thermal,soil,bcu,bch,fluxes]             ...
    = ebal_bigleaf(constants,options,rad,gap,  ...
    meteo,soil,canopy,leafbio)
% function ebal.m calculates the energy balance of a vegetated surface
%
% authors:      Christiaan van der Tol (c.vandertol@utwente.nl)
%               Joris Timmermans 
% date          26 Nov 2007 (CvdT)
% updates       29 Jan 2008 (JT & CvdT)     converted into a function
%               11 Feb 2008 (JT & CvdT)     improved soil heat flux and temperature calculation
%               14 Feb 2008 (JT)            changed h in to hc (as h=Avogadro`s constant)
%               31 Jul 2008 (CvdT)          Included Pntot in output
%               19 Sep 2008 (CvdT)          Converted F0 and F1 from units per aPAR into units per iPAR
%               07 Nov 2008 (CvdT)          Changed layout
%               18 Sep 2012 (CvdT)          Changed Oc, Cc, ec
%                  Feb 2012 (WV)            introduced structures for variables
%                  Sep 2013 (JV, CvT)       introduced additional biochemical model
%               10 Dec 2019 (CvdT)          made a light version (layer
%                                           averaged fluxes)
% parent: SCOPE.m (script)
% uses:
%       RTMt_sb.m, RTMt_planck.m (optional), RTMf.m (optional)
%       resistances.m
%       heatfluxes.m
%       biochemical.m
%       soil_respiration.m
%
% Table of contents of the function
%
%   1. Initialisations for the iteration loop
%           intial values are attributed to variables
%   2. Energy balance iteration loop
%           iteration between thermal RTM and surface fluxes
%   3. Write warnings whenever the energy balance did not close
%   4. Calculate vertical profiles (optional)
%   5. Calculate spectrally integrated energy, water and CO2 fluxes
%
% The energy balance iteration loop works as follows:
%
% RTMo              More or less the classic SAIL model for Radiative
%                   Transfer of sun and sky light (no emission by the vegetation)
% While continue	Here an iteration loop starts to close the energy
%                   balance, i.e. to match the micro-meteorological model
%                   and the radiative transfer model
% 	RTMt_sb         A numerical Radiative Transfer Model for thermal
%                   radiation emitted by the vegetation
% 	resistances     Calculates aerodynamic and boundary layer resistances
%                   of vegetation and soil (the micro-meteorological model)
% 	biochemical     Calculates photosynthesis, fluorescence and stomatal
%                   resistance of leaves (or biochemical_MD12: alternative)
% 	heatfluxes      Calculates sensible and latent heat flux of soil and
%                   vegetation
%                   Next soil heat flux is calculated, the energy balance
%                   is evaluated, and soil and leaf temperatures adjusted
%                   to force energy balance closure
% end {while continue}
%
% meanleaf          Integrates the fluxes over all leaf inclinations
%                   azimuth angles and layers, integrates over the spectrum
%
% usage:
%[iter,fluxes,rad,profiles,thermal]             ...
%         = ebal(iter,options,spectral,rad,gap,leafopt,  ...
%                angles,meteo,soil,canopy,leafbio)
%
% The input and output are structures. These structures are further
% specified in a readme file.
%
% Input:
%
%   iter        numerical parameters used in the iteration for energy balance closure
%   options     calculation options
%   spectral    spectral resolutions and wavelengths
%   rad         incident radiation
%   gap         probabilities of direct light penetration and viewing
%   leafopt     leaf optical properties
%   angles      viewing and observation angles
%   soil        soil properties
%   canopy      canopy properties
%   leafbio     leaf biochemical parameters
%
% Output:
%
%   iter        numerical parameters used in the iteration for energy balance closure
%   fluxes      energy balance, turbulent, and CO2 fluxes
%   rad         radiation spectra
%   thermal     temperatures, aerodynamic resistances and friction velocity
%   bcu, bch    leaf biochemical outputs for sunlit and shaded leaves,
%               respectively

%% 1. initialisations and other preparations for the iteration loop
% initialisations
counter     = 0;              %           Iteration counter of ebal
maxit       = 400;
maxEBer     = 1;
Wc          = 1;
CONT        = 1;              %           is 0 when the calculation has finished

% functions for saturated vapour pressure
es_fun      = @(T)6.107*10.^(7.5.*T./(237.3+T));
s_fun       = @(es, T) es*2.3026*7.5*237.3./(237.3+T).^2;

Ta          = meteo.Ta;
ea          = meteo.ea;
Ca          = meteo.Ca;
Ts          = ones(2,1)*meteo.Ta;
p           = meteo.p;

nl = canopy.nlayers;

Rnuc        = rad.Rnuc;
Tch         = (Ta+.1)*ones(nl,1);                 % Leaf temperature (shaded leaves)
Tcu         = (Ta+.3)*ones(size(Rnuc));           % Leaf tempeFrature (sunlit leaves)
ech         = ea*ones(nl,1);          % Leaf boundary vapour pressure (shaded/sunlit leaves)
Cch         = Ca*ones(nl,1);
ecu         = ea+0*Rnuc;
Ccu         = Ca+0*Rnuc;          % Leaf boundary CO2 (shaded/sunlit leaves)
meteo.L     = -1;                           % Monin-Obukhov length

MH2O        = constants.MH2O;
Mair        = constants.Mair;
rhoa        = constants.rhoa;
cp          = constants.cp;
sigmaSB     = constants.sigmaSB;
Ps          = gap.Ps;
nl          = canopy.nlayers;
kV          = canopy.kV;
xl          = canopy.xl;
LAI         = canopy.LAI;

% other preparations
e_to_q      = MH2O/Mair./p;             % Conversion of vapour pressure [Pa] to absolute humidity [kg kg-1]
Fs          = [1-Ps(end),Ps(end)];      % Matrix containing values for 1-Ps and Ps of soil
Fc          = sum(1-Ps(1:end-1))'/nl;      % Matrix containing values for Ps of canopy
fV          = exp(kV*xl(1:end-1));       % Vertical profile of Vcmax

if size(Rnuc,2)>1
    fVu      = ones(13,36,nl);
    for i = 1:nl
        fVu(:,:,i) = fV(i);
    end
else
    fVu = fV;
end

%% 2. Energy balance iteration loop

%Energy balance loop (Energy balance and radiative transfer)

while CONT                          % while energy balance does not close
    % 2.1. Net radiation of the components
    % Thermal radiative transfer model for vegetation emission (with Stefan-Boltzman's equation)
    rad  = RTMt_sb(constants,rad,soil,leafbio,canopy,gap,Tcu,Tch,Ts(2),Ts(1),0);
    Rnhc    = rad.Rnhc + rad.Rnhct;             %           Canopy (shaded) net radiation
    Rnuc    = rad.Rnuc + rad.Rnuct;             %           Canopy (sunlit) net radiation
    Rnhs    = rad.Rnhs+rad.Rnhst;             %           Soil   (sun+sh) net radiation
    Rnus    = rad.Rnus+rad.Rnust;
    Rns     = [Rnhs Rnus]';
    
 %    rad2  = RTMt_sb(constants,rad,soil,leafbio,canopy,gap,Tcu,Tch,Ts(2),Ts(1),0);
    
    % 2.2. Aerodynamic roughness
    % calculate friction velocity [m s-1] and aerodynamic resistances [s m-1]
    [resist_out]  = resistances(constants,soil,canopy,meteo);
    
    meteo.ustar = resist_out.ustar;
    raa     = resist_out.raa;
    rawc    = resist_out.rawc;
    raws    = resist_out.raws;
    
    % 2.3. Biochemical processes
    [meteo_h,meteo_u]  = deal(meteo);
    meteo_h.T       = mean(Tch)*(1-Fc) + mean(Tcu(:))*Fc;
    meteo_h.eb      = mean(ech)*(1-Fc) + mean(ecu(:))*Fc;
    meteo_h.Cs      = mean(Cch)*(1-Fc) + mean(Cch(:))*Fc;
    meteo_h.Q       = mean(rad.Pnh_Cab)*(1-Fc) + mean(rad.Pnu_Cab(:))*Fc;
%     meteo_u.T       = mean(Tcu(:));
%     meteo_u.eb      = mean(ecu(:));
%     meteo_u.Cs      = mean(Ccu(:));
%     meteo_u.Q       = mean(rad.Pnu_Cab(:));
    if options.Fluorescence_model == 2
        b       = @biochemical_MD12;
    else 
        b       = @biochemical;
    end
%     bch     = b(leafbio,meteo_h,options,constants,fV);
%     bcu     = b(leafbio,meteo_u,options,constants,fVu);
    bch     = b(leafbio,meteo_h,options,constants,1);
%     bcu     = b(leafbio,meteo_u,options,constants,1);
    
    % 2.4. Fluxes (latent heat flux (lE), sensible heat flux (H) and soil heat flux G
    % in analogy to Ohm's law, for canopy (c) and soil (s). All in units of [W m-2]
    
    rss     = soil.rss;
    rac     = (LAI+1)*(raa+rawc);
    ras     = (LAI+1)*(raa+raws);
    [lEch,Hch,ech,Cch,lambdah,sh]     = heatfluxes(rac,bch.rcw,meteo_h.T,ea,Ta,e_to_q,Ca,bch.Ci,constants, es_fun, s_fun);
%     [lEcu,Hcu,ecu,Ccu,lambdau,su]     = heatfluxes(rac,bcu.rcw,meteo_u.T,ea,Ta,e_to_q,Ca,bcu.Ci,constants, es_fun, s_fun);
    [lEs,Hs,~,~,lambdas,ss]           = heatfluxes(ras,rss ,Ts ,ea,Ta,e_to_q,Ca,Ca,constants, es_fun, s_fun);
    
    % integration over the layers and sunlit and shaded fractions
    Hstot   = Fs*Hs;
%     if size(Hcu,2)>1
%             Hctot   = LAI*(Fc*Hch + meanleaf(canopy,Hcu,'angles_and_layers',Ps));
%     else
%         Hctot       = LAI*(Fc*Hch + (1-Fc)*Hcu);
%     end
    Hctot = LAI*(Fc*Hch);
    Htot    = Hstot + Hctot;
    
    % ground heat flux
    G       = 0.35*Rns;
    
    % 2.5. Monin-Obukhov length L
    meteo.L     = Monin_Obukhov(constants,meteo,Htot);                                          % [1]
    
    % 2.6. energy balance errors, continue criterion and iteration counter
    EBerch  = mean(Rnhc) -lEch -Hch;
%     EBercu  = mean(Rnuc(:)) -lEcu -Hcu;
    EBers   = Rns  -lEs  -Hs - G;
    
    counter     = counter+1;                   %        Number of iterations
%     maxEBercu   = max(max(max(abs(EBercu))));
    maxEBerch   = max(abs(EBerch));
    maxEBers    = max(abs(EBers));
    
    CONT        = (... % maxEBercu >   maxEBer    |...
        maxEBerch >   maxEBer     |...
        maxEBers  >   maxEBer)    &...
        counter   <   maxit+1;%        Continue iteration?
    if ~CONT
        if any(isnan([maxEBerch, maxEBers]))
            fprintf('WARNING: NaN in fluxes, counter = %i\n', counter)
        end
        break
    end
    if counter==10, Wc = 0.8;  end
    if counter==60; Wc = 0.6;  end
    if counter==100; Wc = 0.2;  end
    
    % 2.7. New estimates of soil (s) and leaf (c) temperatures, shaded (h) and sunlit (1)
    Tch         = Tch + Wc*EBerch./((rhoa*cp)./rac + rhoa*lambdah*e_to_q.*sh./(rac+bch.rcw)+ 4*leafbio.emis*sigmaSB*(meteo_h.T+273.15).^3);
    Tcu         = mean(Tch) .* ones(size(Rnuc)); %+ Wc*EBercu./((rhoa*cp)./rac + rhoa*lambdau*e_to_q.*su./(rac+bcu.rcw)+ 4*leafbio.emis*sigmaSB*(meteo_u.T+273.15).^3);
    Ts          = Ts + Wc*EBers./(rhoa*cp./ras + rhoa*lambdas*e_to_q.*ss/(ras+rss)+ 4*(1-soil.rs_thermal)*sigmaSB*(Ts+273.15).^3);  
    Tch(abs(Tch)>100) = Ta;
%     Tcu(abs(Tcu)>100) = Ta;
end

%% 3. Print warnings whenever the energy balance could not be solved
if counter>=maxit
    fprintf('WARNING: maximum number of iteratations exceeded\n');
    fprintf('Energy balance error sunlit vegetation = %4.2f W m-2\n' ,maxEBercu);
    fprintf('Energy balance error shaded vegetation = %4.2f W m-2\n' ,maxEBerch);
    fprintf('Energy balance error soil              = %4.2f W m-2\n' ,maxEBers);
end

%% 4. some more outputs
iter.counter    = counter;
% this is not yet written to output but making it avg here breaks RTMt_sb
% thermal.Tcu     = mean(Tcu(:));
% thermal.Tch     = mean(Tch);
thermal.Tcu     = mean(Tch) .* ones(size(Rnuc));
thermal.Tch     = Tch;
thermal.Tsu     = Ts(2);
thermal.Tsh     = Ts(1);

fluxes.Rnctot = mean(Rnhc)*(1-Fc) + mean(Rnuc(:))*Fc; %LAI* aggregator(mean(Rnhc), mean(Rnuc(:)), Fc, Ps, canopy);     % net radiation leaves
fluxes.lEctot = LAI* (Fc*lEch); %aggregator(lEch, lEcu, Fc, Ps, canopy);     % latent heat leaves
fluxes.Hctot  = LAI* (Fc*Hch); %aggregator(Hch, Hcu, Fc, Ps, canopy);       % sensible heat leaves
fluxes.Actot  = LAI* (Fc*bch.A); %aggregator(bch.A, bcu.A, Fc, Ps, canopy);   % photosynthesis leaves
fluxes.Tcave  = aggregator(mean(Tch), mean(Tcu(:)), Fc, Ps, canopy);            % mean leaf temperature

fluxes.Rnstot = Fs*Rns;           % Net radiation soil
fluxes.lEstot = Fs*lEs;           % Latent heat soil
fluxes.Hstot  = Fs*Hs;            % Sensible heat soil
fluxes.Gtot   = Fs*G;             % Soil heat flux
fluxes.Tsave  = Fs*Ts;            % Soil temperature
% fluxes.Resp   = Fs*equations.soil_respiration(Ts); %  Soil respiration = 0

fluxes.Rntot = fluxes.Rnctot + fluxes.Rnstot;
fluxes.lEtot = fluxes.lEctot + fluxes.lEstot;
fluxes.Htot = fluxes.Hctot + fluxes.Hstot;
fluxes.rss = rss;

%% faking back layer structure of SCOPE

bch     = b(leafbio,meteo_h,options,constants,fV);
bcu     = b(leafbio,meteo_h,options,constants,fVu);  % notice meteo_h, not _u

end


function flux_tot = aggregator(shaded_flux, sunlit_flux, Fc, Ps, canopy)
    if size(sunlit_flux, 2) > 1
        flux_tot = Fc*shaded_flux + meanleaf(canopy,sunlit_flux,'angles_and_layers',Ps);
    else
        flux_tot = Fc*shaded_flux + (1-Fc)*sunlit_flux;
    end
end

// src/fluxes/ebal_sunshade.m

function [iter,rad,thermal,soil,bcu,bch,fluxes]             ...
    = ebal_sunshade(constants,options,rad,gap,  ...
    meteo,soil,canopy,leafbio)
% function ebal.m calculates the energy balance of a vegetated surface
%
% authors:      Christiaan van der Tol (c.vandertol@utwente.nl)
%               Joris Timmermans 
% date          26 Nov 2007 (CvdT)
% updates       29 Jan 2008 (JT & CvdT)     converted into a function
%               11 Feb 2008 (JT & CvdT)     improved soil heat flux and temperature calculation
%               14 Feb 2008 (JT)            changed h in to hc (as h=Avogadro`s constant)
%               31 Jul 2008 (CvdT)          Included Pntot in output
%               19 Sep 2008 (CvdT)          Converted F0 and F1 from units per aPAR into units per iPAR
%               07 Nov 2008 (CvdT)          Changed layout
%               18 Sep 2012 (CvdT)          Changed Oc, Cc, ec
%                  Feb 2012 (WV)            introduced structures for variables
%                  Sep 2013 (JV, CvT)       introduced additional biochemical model
%               10 Dec 2019 (CvdT)          made a light version (layer
%                                           averaged fluxes)
% parent: SCOPE.m (script)
% uses:
%       RTMt_sb.m, RTMt_planck.m (optional), RTMf.m (optional)
%       resistances.m
%       heatfluxes.m
%       biochemical.m
%       soil_respiration.m
%
% Table of contents of the function
%
%   1. Initialisations for the iteration loop
%           intial values are attributed to variables
%   2. Energy balance iteration loop
%           iteration between thermal RTM and surface fluxes
%   3. Write warnings whenever the energy balance did not close
%   4. Calculate vertical profiles (optional)
%   5. Calculate spectrally integrated energy, water and CO2 fluxes
%
% The energy balance iteration loop works as follows:
%
% RTMo              More or less the classic SAIL model for Radiative
%                   Transfer of sun and sky light (no emission by the vegetation)
% While continue	Here an iteration loop starts to close the energy
%                   balance, i.e. to match the micro-meteorological model
%                   and the radiative transfer model
% 	RTMt_sb         A numerical Radiative Transfer Model for thermal
%                   radiation emitted by the vegetation
% 	resistances     Calculates aerodynamic and boundary layer resistances
%                   of vegetation and soil (the micro-meteorological model)
% 	biochemical     Calculates photosynthesis, fluorescence and stomatal
%                   resistance of leaves (or biochemical_MD12: alternative)
% 	heatfluxes      Calculates sensible and latent heat flux of soil and
%                   vegetation
%                   Next soil heat flux is calculated, the energy balance
%                   is evaluated, and soil and leaf temperatures adjusted
%                   to force energy balance closure
% end {while continue}
%
% meanleaf          Integrates the fluxes over all leaf inclinations
%                   azimuth angles and layers, integrates over the spectrum
%
% usage:
%[iter,fluxes,rad,profiles,thermal]             ...
%         = ebal(iter,options,spectral,rad,gap,leafopt,  ...
%                angles,meteo,soil,canopy,leafbio)
%
% The input and output are structures. These structures are further
% specified in a readme file.
%
% Input:
%
%   iter        numerical parameters used in the iteration for energy balance closure
%   options     calculation options
%   spectral    spectral resolutions and wavelengths
%   rad         incident radiation
%   gap         probabilities of direct light penetration and viewing
%   leafopt     leaf optical properties
%   angles      viewing and observation angles
%   soil        soil properties
%   canopy      canopy properties
%   leafbio     leaf biochemical parameters
%
% Output:
%
%   iter        numerical parameters used in the iteration for energy balance closure
%   fluxes      energy balance, turbulent, and CO2 fluxes
%   rad         radiation spectra
%   thermal     temperatures, aerodynamic resistances and friction velocity
%   bcu, bch    leaf biochemical outputs for sunlit and shaded leaves,
%               respectively

%% 1. initialisations and other preparations for the iteration loop
% initialisations
counter     = 0;              %           Iteration counter of ebal
maxit       = 400;
maxEBer     = 1;
Wc          = 1;
CONT        = 1;              %           is 0 when the calculation has finished

% functions for saturated vapour pressure
es_fun      = @(T)6.107*10.^(7.5.*T./(237.3+T));
s_fun       = @(es, T) es*2.3026*7.5*237.3./(237.3+T).^2;

Ta          = meteo.Ta;
ea          = meteo.ea;
Ca          = meteo.Ca;
Ts          = ones(2,1)*meteo.Ta;
p           = meteo.p;

nl = canopy.nlayers;

Rnuc        = rad.Rnuc;
Tch         = (Ta+.1)*ones(nl,1);                 % Leaf temperature (shaded leaves)
Tcu         = (Ta+.3)*ones(size(Rnuc));           % Leaf tempeFrature (sunlit leaves)
ech         = ea*ones(nl,1);          % Leaf boundary vapour pressure (shaded/sunlit leaves)
Cch         = Ca*ones(nl,1);
ecu         = ea+0*Rnuc;
Ccu         = Ca+0*Rnuc;          % Leaf boundary CO2 (shaded/sunlit leaves)
meteo.L     = -1;                           % Monin-Obukhov length

MH2O        = constants.MH2O;
Mair        = constants.Mair;
rhoa        = constants.rhoa;
cp          = constants.cp;
sigmaSB     = constants.sigmaSB;
Ps          = gap.Ps;
nl          = canopy.nlayers;
kV          = canopy.kV;
xl          = canopy.xl;
LAI         = canopy.LAI;

% other preparations
e_to_q      = MH2O/Mair./p;             % Conversion of vapour pressure [Pa] to absolute humidity [kg kg-1]
Fs          = [1-Ps(end),Ps(end)];      % Matrix containing values for 1-Ps and Ps of soil
Fc          = sum(1-Ps(1:end-1))'/nl;      % Matrix containing values for Ps of canopy
fV          = exp(kV*xl(1:end-1));       % Vertical profile of Vcmax

if size(Rnuc,2)>1
    fVu      = ones(13,36,nl);
    for i = 1:nl
        fVu(:,:,i) = fV(i);
    end
else
    fVu = fV;
end

%% 2. Energy balance iteration loop

%Energy balance loop (Energy balance and radiative transfer)

while CONT                          % while energy balance does not close
    % 2.1. Net radiation of the components
    % Thermal radiative transfer model for vegetation emission (with Stefan-Boltzman's equation)
    rad  = RTMt_sb(constants,rad,soil,leafbio,canopy,gap,Tcu,Tch,Ts(2),Ts(1),0);
    Rnhc    = rad.Rnhc + rad.Rnhct;             %           Canopy (shaded) net radiation
    Rnuc    = rad.Rnuc + rad.Rnuct;             %           Canopy (sunlit) net radiation
    Rnhs    = rad.Rnhs+rad.Rnhst;             %           Soil   (sun+sh) net radiation
    Rnus    = rad.Rnus+rad.Rnust;
    Rns     = [Rnhs Rnus]';
    
 %    rad2  = RTMt_sb(constants,rad,soil,leafbio,canopy,gap,Tcu,Tch,Ts(2),Ts(1),0);
    
    % 2.2. Aerodynamic roughness
    % calculate friction velocity [m s-1] and aerodynamic resistances [s m-1]
    [resist_out]  = resistances(constants,soil,canopy,meteo);
    
    meteo.ustar = resist_out.ustar;
    raa     = resist_out.raa;
    rawc    = resist_out.rawc;
    raws    = resist_out.raws;
    
    % 2.3. Biochemical processes
    [meteo_h,meteo_u]  = deal(meteo);
    meteo_h.T       = mean(Tch);
    meteo_h.eb      = mean(ech);
    meteo_h.Cs      = mean(Cch);
    meteo_h.Q       = mean(rad.Pnh_Cab);
    meteo_u.T       = mean(Tcu(:));
    meteo_u.eb      = mean(ecu(:));
    meteo_u.Cs      = mean(Ccu(:));
    meteo_u.Q       = mean(rad.Pnu_Cab(:));
    if options.Fluorescence_model == 2
        b       = @biochemical_MD12;
    else 
        b       = @biochemical;
    end
%     bch     = b(leafbio,meteo_h,options,constants,fV);
%     bcu     = b(leafbio,meteo_u,options,constants,fVu);
    bch     = b(leafbio,meteo_h,options,constants,1);
    bcu     = b(leafbio,meteo_u,options,constants,1);
    
    % 2.4. Fluxes (latent heat flux (lE), sensible heat flux (H) and soil heat flux G
    % in analogy to Ohm's law, for canopy (c) and soil (s). All in units of [W m-2]
    
    rss     = soil.rss;
    rac     = (LAI+1)*(raa+rawc);
    ras     = (LAI+1)*(raa+raws);
    [lEch,Hch,ech,Cch,lambdah,sh]     = heatfluxes(rac,bch.rcw,meteo_h.T,ea,Ta,e_to_q,Ca,bch.Ci,constants, es_fun, s_fun);
    [lEcu,Hcu,ecu,Ccu,lambdau,su]     = heatfluxes(rac,bcu.rcw,meteo_u.T,ea,Ta,e_to_q,Ca,bcu.Ci,constants, es_fun, s_fun);
    [lEs,Hs,~,~,lambdas,ss]           = heatfluxes(ras,rss ,Ts ,ea,Ta,e_to_q,Ca,Ca,constants, es_fun, s_fun);
    
    % integration over the layers and sunlit and shaded fractions
    Hstot   = Fs*Hs;
    if size(Hcu,2)>1
            Hctot   = LAI*(Fc*Hch + meanleaf(canopy,Hcu,'angles_and_layers',Ps));
    else
        Hctot       = LAI*(Fc*Hch + (1-Fc)*Hcu);
    end
    Htot    = Hstot + Hctot;
    
    % ground heat flux
    G       = 0.35*Rns;
    
    % 2.5. Monin-Obukhov length L
    meteo.L     = Monin_Obukhov(constants,meteo,Htot);                                          % [1]
    
    % 2.6. energy balance errors, continue criterion and iteration counter
    EBerch  = mean(Rnhc) -lEch -Hch;
    EBercu  = mean(Rnuc(:)) -lEcu -Hcu;
    EBers   = Rns  -lEs  -Hs - G;
    
    counter     = counter+1;                   %        Number of iterations
    maxEBercu   = max(max(max(abs(EBercu))));
    maxEBerch   = max(abs(EBerch));
    maxEBers    = max(abs(EBers));
    
    CONT        = ( maxEBercu >   maxEBer    |...
        maxEBerch >   maxEBer     |...
        maxEBers  >   maxEBer)    &...
        counter   <   maxit+1;%        Continue iteration?
    if ~CONT
        if any(isnan([maxEBercu, maxEBerch, maxEBers]))
            fprintf('WARNING: NaN in fluxes, counter = %i\n', counter)
        end
        break
    end
    if counter==10, Wc = 0.8;  end
    if counter==60; Wc = 0.6;  end
    if counter==100; Wc = 0.2;  end
    
    % 2.7. New estimates of soil (s) and leaf (c) temperatures, shaded (h) and sunlit (1)
    Tch         = Tch + Wc*EBerch./((rhoa*cp)./rac + rhoa*lambdah*e_to_q.*sh./(rac+bch.rcw)+ 4*leafbio.emis*sigmaSB*(meteo_h.T+273.15).^3);
    Tcu         = Tcu + Wc*EBercu./((rhoa*cp)./rac + rhoa*lambdau*e_to_q.*su./(rac+bcu.rcw)+ 4*leafbio.emis*sigmaSB*(meteo_u.T+273.15).^3);
    Ts          = Ts + Wc*EBers./(rhoa*cp./ras + rhoa*lambdas*e_to_q.*ss/(ras+rss)+ 4*(1-soil.rs_thermal)*sigmaSB*(Ts+273.15).^3);  
    Tch(abs(Tch)>100) = Ta;
    Tcu(abs(Tcu)>100) = Ta;
end

%% 3. Print warnings whenever the energy balance could not be solved
if counter>=maxit
    fprintf('WARNING: maximum number of iteratations exceeded\n');
    fprintf('Energy balance error sunlit vegetation = %4.2f W m-2\n' ,maxEBercu);
    fprintf('Energy balance error shaded vegetation = %4.2f W m-2\n' ,maxEBerch);
    fprintf('Energy balance error soil              = %4.2f W m-2\n' ,maxEBers);
end

%% 4. some more outputs
iter.counter    = counter;
% this is not yet written to output but making it avg here breaks RTMt_sb
% thermal.Tcu     = mean(Tcu(:));
% thermal.Tch     = mean(Tch);
thermal.Tcu     = Tcu;
thermal.Tch     = Tch;
thermal.Tsu     = Ts(2);
thermal.Tsh     = Ts(1);

fluxes.Rnctot = LAI* aggregator(mean(Rnhc), mean(Rnuc(:)), Fc, Ps, canopy);     % net radiation leaves
fluxes.lEctot = LAI* aggregator(lEch, lEcu, Fc, Ps, canopy);     % latent heat leaves
fluxes.Hctot  = LAI* aggregator(Hch, Hcu, Fc, Ps, canopy);       % sensible heat leaves
fluxes.Actot  = LAI* aggregator(bch.A, bcu.A, Fc, Ps, canopy);   % photosynthesis leaves
fluxes.Tcave  = aggregator(mean(Tch), mean(Tcu(:)), Fc, Ps, canopy);            % mean leaf temperature

fluxes.Rnstot = Fs*Rns;           % Net radiation soil
fluxes.lEstot = Fs*lEs;           % Latent heat soil
fluxes.Hstot  = Fs*Hs;            % Sensible heat soil
fluxes.Gtot   = Fs*G;             % Soil heat flux
fluxes.Tsave  = Fs*Ts;            % Soil temperature
% fluxes.Resp   = Fs*equations.soil_respiration(Ts); %  Soil respiration = 0

fluxes.Rntot = fluxes.Rnctot + fluxes.Rnstot;
fluxes.lEtot = fluxes.lEctot + fluxes.lEstot;
fluxes.Htot = fluxes.Hctot + fluxes.Hstot;
fluxes.rss = rss;

%% faking back layer structure of SCOPE

bch     = b(leafbio,meteo_h,options,constants,fV);
bcu     = b(leafbio,meteo_u,options,constants,fVu);

end


function flux_tot = aggregator(shaded_flux, sunlit_flux, Fc, Ps, canopy)
    if size(sunlit_flux, 2) > 1
        flux_tot = Fc*shaded_flux + meanleaf(canopy,sunlit_flux,'angles_and_layers',Ps);
    else
        flux_tot = Fc*shaded_flux + (1-Fc)*sunlit_flux;
    end
end

// src/fluxes/heatfluxes.m

function [lE, H, ec, Cc, lambda, s]  = heatfluxes(ra,rs,Tc,ea,Ta,e_to_q,Ca,Ci,constants, es_fun, s_fun)    

% author: Dr. ir. Christiaan van der Tol (c.vandertol@utwente.nl)
% date:     7 Dec 2007
% updated: 15 Apr 2009 CvdT     changed layout
% updated: 14 Sep 2012 CvdT     added ec and Cc to output
% updated: 09 Dec 2019 CvdT     modified for computational efficiency
%
% parent: ebal.m
%
% usage:
% function [lE, H]  = heatfluxes(ra,rs,Tc,ea,Ta,e_to_q,PSI,Ca,Ci,constants,es_fun, s_fun)
% 
% this function calculates latent and sensible heat flux
%
% input:
%   ra          aerodynamic resistance for heat         s m-1
%   rs          stomatal resistance                     s m-1
%   Tc          leaf temperature                        oC
%   ea          vapour pressure above canopy            hPa
%   Ta          air temperature above canopy            oC
%   e_to_q      conv. from vapour pressure to abs hum   hPa-1
%   PSI         leaf water potential                    J kg-1
%   Ca          ambient CO2 concentration               umol m-3
%   Ci          intercellular CO2 concentration         umol m-3
%   constants   a structure with physical constants
%   es_fun      saturated pressure function es(hPa)=f(T(C))
%   s_fun       slope of the saturated pressure function (s(hPa/C) = f(T(C), es(hPa))
%
% output:
%   lEc         latent heat flux of a leaf              W m-2
%   Hc          sensible heat flux of a leaf            W m-2
%   ec          vapour pressure at the leaf surface     hPa
%   Cc          CO2 concentration at the leaf surface   umol m-3

rhoa = constants.rhoa;
cp   = constants.cp;
%MH2O = constants.MH2O;
%R    = constants.R;

lambda      = (2.501-0.002361*Tc)*1E6;  %      [J kg-1]  Evapor. heat (J kg-1)
ei = es_fun(Tc);
s = s_fun(ei, Tc);

%ei          = es.*exp(1E-3*PSI*MH2O/R./(Tc+273.15));  
qi          = ei .* e_to_q;
qa          = ea .* e_to_q;

lE          = rhoa./(ra+rs).*lambda.*(qi-qa);   % [W m-2]   Latent heat flux
H           = (rhoa*cp)./ra.*(Tc-Ta);           % [W m-2]   Sensible heat flux
ec          = ea + (ei-ea).*ra./(ra+rs);         % [W m-2] vapour pressure at the leaf surface
Cc          = Ca - (Ca-Ci).*ra./(ra+rs);        % [umol m-2 s-1] CO2 concentration at the leaf surface
end

// src/fluxes/resistances.m

function [resist_out] = resistances(constants,soil,canopy,meteo)
%
%   function resistances calculates aerodynamic and boundary resistances
%   for soil and vegetation
%
%   Date:       01 Feb 2008
%   Authors:    Anne Verhoef            (a.verhoef@reading.ac.uk)
%               Christiaan van der Tol  (tol@itc.nl)
%               Joris Timmermans        (j_timmermans@itc.nl)
%   Source:     Wallace and Verhoef (2000) 'Modelling interactions in
%               mixed-plant communities: light, water and carbon dioxide', in: Bruce
%               Marshall, Jeremy A. Roberts (ed), 'Leaf Development and Canopy Growth',
%               Sheffield Academic Press, UK. ISBN 0849397693
%               
%               ustar:  Tennekes, H. (1973) 'The logaritmic wind profile', J.
%               Atmospheric Science, 30, 234-238
%               Psih:   Paulson, C.A. (1970), The mathematical
%               representation of wind speed and temperature in the
%               unstable atmospheric surface layer. J. Applied Meteorol. 9,
%               857-861
%       
% Note: Equation numbers refer to equation numbers in Wallace and Verhoef (2000)
% 
% Usage:
%   [resist_out] = resistances(resist_in)
%
% The input and output are structures. These structures are further
% specified in a readme file.
%
% Input: 
%   resist_in   aerodynamic resistance parameters and wind speed
%
% The strucutre resist_in contains the following elements:
% u         =   windspeed
% L         =   stability
% LAI       =   Leaf Area Index

% rbs       =   Boundary Resistance of soil                         [s m-1]
% rss       =   Surface resistance of soil for vapour transport     [s m-1]
% rwc       =   Within canopy Aerodynamic Resistance canopy         [s m-1]

% z0m       =   Roughness lenght for momentum for the vegetation    [m]
% d         =   Displacement height (Zero plane)                    [m]
% z         =   Measurement height                                  [m]
% h         =   Vegetation height                                   [m]

%
% Output:
%   resist_out  aeorodynamic resistances
%
% The strucutre resist_out contains the following elements:
% ustar     =   Friction velocity                                   [m s-1]
% raa       =   Aerodynamic resistance above the canopy             [s m-1]                     
% rawc      =   Total resistance within the canopy (canopy)         [s m-1]
% raws      =   Total resistance within the canopy (soil)           [s m-1]

% rai       =   Aerodynamic resistance in inertial sublayer         [s m-1]
% rar       =   Aerodynamic resistance in roughness sublayer        [s m-1]
% rac       =   Aerodynamic resistance in canopy layer (above z0+d) [s m-1]

% rbc       =   Boundary layer resistance (canopy)                  [s m-1]
% rwc       =   Aerodynamic Resistance within canopy(canopy)(Update)[s m-1]

% rbs       =   Boundary layer resistance (soil) (Update)           [s m-1]
% rws       =   Aerodynamic resistance within canopy(soil)          [s m-1] 

% rss       =   Surface resistance vapour transport(soil)(Update)   [s m-1]

% uz0       =   windspeed at z0                                     [m s-1]
% Kh        =   Diffusivity for heat                                [m2s-1]

%% parameters
%global constants
kappa   = constants.kappa;
Cd      =  canopy.Cd;
LAI     =  canopy.LAI;
rwc     =  canopy.rwc;
z0m     =  canopy.zo;
d       =  canopy.d;
h       =  canopy.hc;
%w       =  canopy.leafwidth;
z       =  meteo.z;
u       =  max(0.3,meteo.u);
L       =  meteo.L;
rbs     =  soil.rbs;
%rss       =  resist_in.rss;

% derived parameters
%zr: top of roughness sublayer, bottom of intertial sublayer
zr			= 2.5*h;                   %                            [m]			
%n: dimensionless wind extinction coefficient                       W&V Eq 33
n			= Cd*LAI/(2*kappa^2);      %                            [] 

%% stability correction for non-neutral conditions
unst        = (L <  0 & L>-500);
st          = (L >  0 & L<500);  
x       	= (1-16*z./L).^(1/4); % only used for unstable

% stability correction functions, friction velocity and Kh=Km=Kv
pm_z    	= psim(z -d,L,unst,st,x);
ph_z    	= psih(z -d,L,unst,st,x);
pm_h        = psim(h -d,L,unst,st,x);
%ph_h       = psih(h -d,L,unst,st);
ph_zr       = psih(zr-d,L,unst,st,x).*(z>=zr) + ph_z.*(z<zr);
phs_zr      = phstar(zr,zr,d,L,st,unst,x);
phs_h		= phstar(h ,zr,d,L,st,unst,x);

ustar   	= max(.001,kappa*u./(log((z-d)/z0m) - pm_z));%          W&V Eq 30
Kh          = kappa*ustar*(zr-d);                  %                W&V Eq 35

if unst
    resist_out.Kh	= Kh*(1-16*(h-d)./L).^.5;% W&V Eq 35
elseif st
    resist_out.Kh   = Kh*(1+ 5*(h-d)./L  ).^-1;% W&V Eq 35
else
    resist_out.Kh = Kh;
end

%% wind speed at height h and z0m
uh			= max(ustar/kappa .* (log((h-d)/z0m) - pm_h     ),.01);
uz0 		= uh*exp(n*((z0m+d)/h-1));                      %       W&V Eq 32

%% resistances

resist_out.uz0 = uz0; 
resist_out.ustar = ustar;
rai = (z>zr).*(1./(kappa*ustar).*(log((z-d) /(zr-d))  - ph_z   + ph_zr));% W&V Eq 41 
rar = 1./(kappa*ustar).*((zr-h)/(zr-d)) 	 - phs_zr + phs_h;% W&V Eq 39
rac = h*sinh(n)./(n*Kh)*(log((exp(n)-1)/(exp(n)+1)) - log((exp(n*(z0m+ d )/h)-1)/(exp(n*(z0m +d )/h)+1))); % W&V Eq 42
rws = h*sinh(n)./(n*Kh)*(log((exp(n*(z0m+d)/h)-1)/(exp(n*(z0m+d)/h)+1)) - log((exp(n*(.01    )/h)-1)/(exp(n*(.01    )/h)+1))); % W&V Eq 43
%rbc = 70/LAI * sqrt(w./uz0);						%		W&V Eq 31, but slightly different

resist_out.rai = rai;
resist_out.rar = rar;
resist_out.rac = rac;
resist_out.rws = rws;
%resist_out.rbc = rbc;

raa  = rai + rar + rac;
rawc = rwc;% + rbc;
raws = rws + rbs;

resist_out.raa  = raa;          % aerodynamic resistance above the canopy           W&V Figure 8.6
resist_out.rawc	= rawc;			% aerodynamic resistance within the canopy (canopy)
resist_out.raws	= raws;			% aerodynamic resistance within the canopy (soil)

return


%% subfunction pm for stability correction (eg. Paulson, 1970)
function pm = psim(z,L,unst,st,x)
pm      	= 0;
if unst
    pm          = 2*log((1+x)/2)+log((1+x.^2)/2) -2*atan(x)+pi/2;   %   unstable
elseif st
    pm          = -5*z./L;                                      %   stable
end
return

%% subfunction ph for stability correction (eg. Paulson, 1970)
function ph = psih(z,L,unst,st,x)
ph        = 0;
if unst
    ph      = 2*log((1+x.^2)/2);                                %   unstable
elseif st
    ph      = -5*z./L;                                      %   stable
end
return

%% subfunction ph for stability correction (eg. Paulson, 1970)
function phs = phstar(z,zR,d,L,st,unst,x)
phs         = 0;
if unst
    phs     = (z-d)/(zR-d)*(x.^2-1)./(x.^2+1);
elseif st
    phs     = -5*z./L;
end
return

// src/plotting/v10_biochemical_one_leaf.m
%% youtube video 10
% https://youtu.be/COM89KOGKMo?si=1xHLTl1sPibPoyCP

%% canopy
fV = 1;  % the exponent of the vertical decline of Vcmax in the canopy

%% options_adj

options_adj = struct();

options_adj.apply_T_corr = 0;

%% constant

constants_adj = struct();

constants_adj.R    = 8.314;      % [J mol-1K-1]  Molar gas constant
constants_adj.rhoa = 1.2047;     % [kg m-3]      Specific mass of air
constants_adj.Mair = 28.96;      % [g mol-1]     Molecular mass of dry air


%% meteo

meteo_leaf = struct();

% 1 W m-2 ~ 4.6 umol photons m-2 s-1
meteo_leaf.Q = 500;  % umol photons m-2 s-1, aPAR by chlorophyll
meteo_leaf.T = 20;     % deg C,     leaf temperature
meteo_leaf.Cs = 410;   % ppm,       leaf boundary layer CO2 concentration
meteo_leaf.eb = 15;    % hPa,       actual atmospheric vapour pressure  meteo.ea
meteo_leaf.Oa = 209;  % per mille, atmospheric O2 concentraion  %meteo.Oa;
meteo_leaf.p = 970;    % hPa,       atmospheric pressure, meteo.p;

%% leaf biochemical parameters

leafbio_adj = struct();

% in Magniani C4 CO2 response of ps, Fs, eta is not present
leafbio_adj.Type = 'C3'; 

leafbio_adj.Vcmax25 = 60;
leafbio_adj.RdPerVcmax25 = 0.015;

leafbio_adj.BallBerrySlope = 8;
leafbio_adj.BallBerry0 = 0.01;


% CvdT model
% Ф(Kp) + Ф(Kf) + Ф(Kd) + Ф(Kn) = 1
% Kn is modelled based on the degree of light saturation: 'x'
leafbio_adj.Kn0 = 2.48;
leafbio_adj.Knalpha = 2.83;
leafbio_adj.Knbeta = 0.1140;

leafbio_adj.stressfactor = 1;  % for Vcmax

leafbio_adj.TDP = define_temp_response_biochem; % temperature response C3 and C4 according to CLM4 model


% Magniani model
% Ф(Kp) + Ф(Kf) + Ф(Kd) + Ф(Kn) = 1
% Kp and Kd are adjusted witht the respective stress factors
leafbio_adj.qLs = 1; % fraction of functional reaction centres (RC)
leafbio_adj.kNPQs = 0; % rate constant of sustained thermal dissipation, normalized to (Kf+Kd)

leafbio_adj.beta = 0.51; % fraction of photons partitioned to PSII

leafbio_adj.Tyear = 15; % mean annual temperature



%% light response curves

meteo_leaf.Q = linspace(1, 2000, 100);
meteo_leaf.Cs = 410;
meteo_leaf.T = 20; 
x_values = meteo_leaf.Q;
x_label = 'aPAR\_Cab, \mumol photons m-2 s-1';

leaf_out_vdt = biochemical(leafbio_adj,meteo_leaf,options_adj,constants_adj,fV);
leaf_out_md12 = biochemical_MD12(leafbio_adj,meteo_leaf,options_adj,constants_adj,fV);

plot_leaf(x_values, x_label, leaf_out_vdt, leaf_out_md12, constants_adj)

%% CO2 response curves

meteo_leaf.Q = 500;
meteo_leaf.Cs = linspace(0, 1000, 100); 
meteo_leaf.T = 20; 
x_values = meteo_leaf.Cs;
x_label = 'Cs, ppm CO_2';

leaf_out_vdt = biochemical(leafbio_adj,meteo_leaf,options_adj,constants_adj,fV);
leaf_out_md12 = biochemical_MD12(leafbio_adj,meteo_leaf,options_adj,constants_adj,fV);

plot_leaf(x_values, x_label, leaf_out_vdt, leaf_out_md12, constants_adj)

%% Temperature response

meteo_leaf.Q = 500;
meteo_leaf.Cs = 410; 
meteo_leaf.T = linspace(-10, 50, 10);
x_values = meteo_leaf.T;
x_label = 'T_{leaf}, ^oC';

leaf_out_vdt = biochemical(leafbio_adj,meteo_leaf,options_adj,constants_adj,fV);
leaf_out_md12 = biochemical_MD12(leafbio_adj,meteo_leaf,options_adj,constants_adj,fV);

plot_leaf(x_values, x_label, leaf_out_vdt, leaf_out_md12, constants_adj)

%% 

function plot_leaf(x_values, x_label, leaf_out_vdt, leaf_out_md12, constants_adj)
    
    figure

    subplot(2,3,1)
    plot(x_values, leaf_out_vdt.A, 'DisplayName', ['Fluorescence_model == 0', ...
        newline, 'Van der Tol et al., 2014'])
    hold on
    plot(x_values, leaf_out_md12.A, 'DisplayName', ['Fluorescence_model == 1', ...
        newline, 'Magniani'])
    
    xlabel(x_label)
    ylabel('A, \mumol CO_2 m-2 s-1')
    title('Net photosynthesis')
    
    % legend('location', 'best', 'Interpreter','none')
    
    
    subplot(2,3,2)
    plot(x_values, leaf_out_vdt.Ja, 'DisplayName', ['Fluorescence_model == 0', ...
        newline, 'Van der Tol et al., 2014'])
    hold on
    plot(x_values, leaf_out_md12.Ja, 'DisplayName', ['Fluorescence_model == 1', ...
        newline, 'Magniani'])
    
    xlabel(x_label)
    ylabel('Ja, \mumol electrons m-2 s-1')
    title('Electron transport rate')
    
    % legend('location', 'best', 'Interpreter','none')
    
    
    
    subplot(2,3,3)
    
    num = constants_adj.rhoa./(constants_adj.Mair*1E-3);
    % gs = (rhoa./(Mair*1E-3)) / rcw  % mol H_2O m-2 s-1
    plot(x_values, num ./ leaf_out_vdt.rcw, 'DisplayName', ['Fluorescence_model == 0', ...
        newline, 'Van der Tol et al., 2014'])
    hold on
    plot(x_values, num./leaf_out_md12.rcw, 'DisplayName', ['Fluorescence_model == 1', ...
        newline, 'Magniani'])
    
    xlabel(x_label)
    % ylabel('rcw, s m-1')
    ylabel('gsw, mol H_2O m-2 s-1')
    title('Stomatal conductance')
    
    % legend('location', 'best', 'Interpreter','none')
    
    
    subplot(2,3,4)
    plot(x_values, leaf_out_vdt.ps, 'DisplayName', ['Fluorescence_model == 0', ...
        newline, 'Van der Tol et al., 2014'])
    hold on
    plot(x_values, leaf_out_md12.ps, 'DisplayName', ['Fluorescence_model == 1', ...
        newline, 'Magniani'])
    
    xlabel(x_label)
    ylabel('ps, \Phi_{PSII}')
    title('PSII photochemcal quantum yield')
    
    % legend('location', 'best', 'Interpreter','none')
    
    
    subplot(2,3,5)
    plot(x_values, leaf_out_vdt.fs, 'DisplayName', ['Fluorescence_model == 0', ...
        newline, 'Van der Tol et al., 2014'])
    hold on
    plot(x_values, leaf_out_md12.fs, 'DisplayName', ['Fluorescence_model == 1', ...
        newline, 'Magniani'])
    
    xlabel(x_label)
    ylabel('fs, \Phi_F_s')
    title('Steady-state fluorescence')
    
    % legend('location', 'best', 'Interpreter','none')
    
    
    subplot(2,3,6)
    plot(x_values, leaf_out_vdt.eta, 'DisplayName', ['Fluorescence_model == 0', ...
        newline, 'Van der Tol et al., 2014'])
    hold on
    plot(x_values, leaf_out_md12.eta, 'DisplayName', ['Fluorescence_model == 1', ...
        newline, 'Magniani'])
    
    xlabel(x_label)
    ylabel('eta, F_s / F_o')
    title('eta, F_s / F_o')
    
    legend('location', 'best', 'Interpreter','none')
    
    set(findall(gcf,'-property','FontSize'),'FontSize', 14)
end

// src/plotting/v5_plot_spectral.m

%% youtube video 5
% https://youtu.be/6DvOiadMA_M?si=OWPUq3KQp2nl0eKi

%% reflectance
% components
% (1) soil
% (2) leaf (per layer)
% (3) canopy

figure

wl = spectral.wlS * 1e-3; 

plot(wl, leafopt.refl(1, :), 'o', 'DisplayName', 'leafopt.refl (leaf)')
hold on
plot(wl, soil.refl, 'o', 'DisplayName', 'soil.refl (soil)')

plot(wl, rad.refl, 'x', 'DisplayName', 'rad.refl (canopy)')

set(findall(gcf,'-property','FontSize'),'FontSize', 14)
set(gca, 'XScale', 'log')
legend

xlabel('wavelength, \mum')
title('reflectance')


%% spectral regions (see MATLAB structure "spectral" in the workspace)
%
% (1) PROSPECT 400 - 2400 nm
% (2) TIR - 2400-50000 nm (spacing 100 nm)
% (3) SIF - 640-850 nm  (also SIF excitation 400-750 nm)
% (4) xanthophyll - 500 - 600 nm

%% spectral regions

figure

w = 10; % line width

plot(spectral.wlS*1e-3, repmat(1,size(spectral.wlS)), 'Color', "#0072BD", 'LineWidth', w)
hold on
plot(spectral.wlP*1e-3, repmat(2,size(spectral.wlP)), 'Color', "#77AC30", 'LineWidth', w)

% plot(spectral.wlE*1e-3, repmat(2.5,size(spectral.wlE)), 'o', 'Color', "#77AC30", 'DisplayName', 'spectral.wlF (SIF excitation)')

plot(spectral.wlF*1e-3, repmat(3,size(spectral.wlF)), 'Color', "#EDB120", 'LineWidth', w)

plot(spectral.wlT*1e-3, repmat(4,size(spectral.wlT)), 'Color', "#D95319", 'LineWidth', w)

plot(spectral.wlZ*1e-3, repmat(5,size(spectral.wlZ)), 'Color', "#7E2F8E", 'LineWidth', w)

ylim([0, 6])
yticks([1,2,3,4,5])
yticklabels({'wlS (SCOPE)','wlP (FLUSPECT)','wlF (SIF)','wlT (TIR)', 'wlZ (PRI)'})

% legend
title('spectral regions of the SCOPE model', 'spectral.wlS')
xlabel('wavelength, \mum')
set(findall(gcf,'-property','FontSize'),'FontSize', 14)
set(gca, 'XScale', 'log')

%% radiance background
% (1) hemispherically integrated flux, E
% (2) flux in the observation direction, L, sr-1
%% L

figure

subplot(1,2,1)

wl = spectral.wlS * 1e-3; 

plot(wl, rad.Lo_, 'DisplayName', ['rad.Lo_' newline 'reflected'])
hold on

plot(wl, rad.Lot_, 'DisplayName', ['rad.Lot_' newline 'emitted in TIR wl'])

wl = spectral.wlF * 1e-3;  % um
plot(wl, rad.LoF_, 'o', 'DisplayName', ['rad.LoF_' newline 'emitted in SIF wl'])


xlabel('wavelength, \mum')
ylabel('L, W m-2 \mum-1 sr-1')
title('Radiance (L)', 'in the observation direction')

legend('Interpreter', 'none')

set(findall(gcf,'-property','FontSize'),'FontSize', 14)
set(gca, 'XScale', 'log')

%% E
% figure
subplot(1,2,2)

wl = spectral.wlS * 1e-3; 

plot(wl, rad.Eout_, 'DisplayName', ['rad.Eout_' newline 'reflected in FLUSPECT wl'])
hold on

plot(wl, rad.Eoutte_, 'DisplayName', ['rad.Eoutte_' newline 'emitted in TIR wl'])

wl = spectral.wlF * 1e-3;  % um
plot(wl, rad.EoutF_, 'o', 'DisplayName', ['rad.EoutF_' newline 'emitted in SIF wl'])


xlabel('wavelength, \mum')
ylabel('E, W m-2 \mum-1')
title('Radiance (E)', 'hemispherically integrated')

legend('Interpreter', 'none')

set(findall(gcf,'-property','FontSize'),'FontSize', 14)
set(gca, 'XScale', 'log')

%% apparent reflectance

figure

subplot(2,2,1)
wl = spectral.wlS * 1e-3; 

plot(wl, rad.refl, 'DisplayName', 'rad.refl (reflectance)')
hold on
plot(wl, rad.reflapp, 'DisplayName', 'rad.reflapp (apparent reflectance)')

set(gca, 'XScale', 'log')

title('canopy reflectances')
xlabel('wavelength, \mum')
legend('Interpreter', 'none')


subplot(2,2,2)
plot(wl, rad.reflapp - rad.refl, 'Color', "#EDB120")

set(gca, 'XScale', 'log')

title('apparent minus true reflectance', 'rad.reflapp - rad.refl', 'Interpreter','none')
xlabel('wavelength, \mum')


subplot(2,2,3)

plot(wl, rad.Lotot_, 'DisplayName', 'rad.Lotot_ (radiance)')
hold on
plot(wl, rad.Lototf_, 'DisplayName', 'rad.Lototf_ (apparent radiance)')

set(gca, 'XScale', 'log')

title('canopy radiance in observation direction')
xlabel('wavelength, \mum')
ylabel('L, W m-2 \mum-1 sr-1')
legend('Interpreter', 'none')

subplot(2,2,4)

plot(wl, rad.Lototf_ - rad.Lotot_, 'Color', "#EDB120")

set(gca, 'XScale', 'log')

title('apparent minus true radiance', 'rad.Lototf_ - rad.Lotot_', 'Interpreter','none')
xlabel('wavelength, \mum')
ylabel('L, W m-2 \mum-1 sr-1')


set(findall(gcf,'-property','FontSize'),'FontSize', 14)


for i=1:4
    subplot(2,2,i)
    xlim([640, 850]*1e-3)
end



%% xanthophyll cycle (leaf reflectance only)
% this contribution is summed with Lo_, so individual plotting is not
% easlity possible

figure

wl = spectral.wlS * 1e-3; 

plot(wl, leafopt.refl(1, :), 'DisplayName', ['leafopt.refl' newline 'leaf reflectance 100% viola'])
hold on
plot(wl, leafopt.reflZ(1, :), 'DisplayName', ['leafopt.reflZ' newline 'leaf reflectance 100% zea [light stress]'])


legend
set(findall(gcf,'-property','FontSize'),'FontSize', 14)
set(gca, 'XScale', 'log')

title('Leaf reflectance', 'violaxanthin cycle')
xlabel('wavelength, \mum')

% xlim([500, 600]*1e-3)


% p = patch([wl; flip(wl)], [leafopt.refl(1, :), fliplr(leafopt.reflZ(1, :))],  'blue', 'DisplayName', 'effect');
% p.FaceColor = '#7E2F8E';
% p.EdgeColor = 'none';

// src/plotting/v6_plot_directional.m

%% youtube video 6
% https://youtu.be/bIuc4kKLAnk?si=ZSJ_37z82vpG4NlG

%% input (to be adjusted by the user)

path_output =  fullfile('..', '..', Output_dir, 'directional');
% path_output =  "../../output/direction_2024-06-10-1150/directional";

lst = dir(path_output);
fnames = {lst.name};

Angles = load(fullfile(path_output, fnames{contains(fnames, 'Angles')}));

% optical
wl_target = 555; 
Values = load(fullfile(path_output, fnames{contains(fnames, 'refl')}));
% 
% % SIF
% wl_target = 760;  
% Values = load(fullfile(path_output, fnames{contains(fnames, 'Fluorescence')}));
% 
% % TIR
% wl_target = 4000;  
% Values = load(fullfile(path_output, fnames{contains(fnames, 'Thermal radiances')}));


%% grid creation
obs_zenith                  = Angles(1,:);
obs_azimuth               	= Angles(2,:);

obs_azimuth180              = obs_azimuth-360*(obs_azimuth>180);

obs_zenith_i                = 0:90;
obs_azimuth_i               = -180:1:180;

[Obs_Zenith_i,Obs_Azimuth_i]= meshgrid(obs_zenith_i,obs_azimuth_i);

% conversion to polar coordinates
% https://nl.mathworks.com/matlabcentral/answers/716343-creating-polar-mesh-in-matlab#answer_597478
% r = obs_zenith * pi/180
% phi = obs_azimuth180  * pi/180 + pi/2

x                           = obs_zenith  *pi/180 .* cos(obs_azimuth180  *pi/180+pi/2);
y                           = obs_zenith  *pi/180 .* sin(obs_azimuth180  *pi/180+pi/2);

X_i                         = Obs_Zenith_i*pi/180 .* cos(Obs_Azimuth_i*pi/180+pi/2);
Y_i                         = Obs_Zenith_i*pi/180 .* sin(Obs_Azimuth_i*pi/180+pi/2);


%% wavelength selection

wl = Values(:,1    );  % nm
i_wl = find(wl == wl_target);

assert(~isempty(i_wl), 'wavelength %d nm was not found. Adjust wl_target?', wl_target)

BRDF                        = Values(:,2:end);

% regridding 

BRDF_i                      = griddata(x,y,BRDF(i_wl,:),X_i,Y_i,'v4');


%% 2d plot

figure

xli = .5*pi*[0 -1.15  -.1 1];
yli = .5*pi*[1  0 -1.05 0];


subplot(2, 2, 1)
z = pcolor(X_i,Y_i,BRDF_i); 
hold on
set(z,'LineStyle','none')
for k = 1:4
    plot(20*k/180*pi.*cos((-1:.01:1)*pi),20*k/180*pi.*sin((-1:.01:1)*pi),'Color', 'black')
    text(20*k/180*pi-.2,.2,num2str(20*k),'FontSize',12.727272727272727,'Color', 'black');
    text(xli(k),yli(k),num2str(90*(k-1)),'FontSize',14,'Color','k','FontAngle','italic');
end

title(sprintf('BRDF @ %.2d nm', wl_target),'FontSize',14)
colorbar
colormap jet

axis off

%% 3d plot


% figure
subplot(2, 2, 2)
surf(X_i,Y_i,BRDF_i); 
shading interp
colormap jet
colorbar

title(sprintf('BRDF @ %.2d nm', wl_target),'FontSize',14)


%% principal plane: psi 0 -> 180


% figure
subplot(2, 2, 3)


% forward scattering
i = obs_azimuth == 180;
x_f = obs_zenith(i);
y_f = BRDF(i_wl,i);

% back scattering
i = obs_azimuth == 0;
x_b = obs_zenith(i) * -1;
y_b = BRDF(i_wl,i);

x = [x_f, x_b];
[~, i2] = sort(x);
y = [y_f, y_b];

plot(x(i2), y(i2))
xlabel('observation zenith angle (tts)')

title(sprintf('BRDF @ %.2d nm\nprincipal plane (phi 0 -> 180)\n"vertical" where the Sun is', wl_target),'FontSize',14)


%% cross-principal plane: psi 90 -> 270


% figure
subplot(2, 2, 4)

% forward scattering
i = obs_azimuth == 270;
x_f = obs_zenith(i);
y_f = BRDF(i_wl,i);

% back scattering
i = obs_azimuth == 90;
x_b = obs_zenith(i) * -1;
y_b = BRDF(i_wl,i);

x = [x_f, x_b];
[~, i2] = sort(x);
y = [y_f, y_b];

plot(x(i2), y(i2))
xlabel('observation zenith angle (tts)')

title(sprintf('BRDF @ %.2d nm\ncross-principal plane (phi 90 -> 270)\n perpendicular to principal', wl_target),'FontSize',14)

// src/plotting/v7_1_plot_profiles_scalars.m

%% youtube video 7
% https://youtu.be/lhmP78f5tYg?si=O2i3m-WipOSlY8xP

%% plot profiles
% this script works with the structures available in the workspace after
% the SCOPE run

%% contribution of sunlit and shaded leaves depends on their fraction
% computed as the weighted [by the gap fraction gap.Ps] sum L386 RTMo 
% profiles.Pn1d  = ((1-Ps(1:nl)).*Pnhc  + Ps(1:nl).*(Pnu1d));  

nl = canopy.nlayers;
Ps = gap.Ps;

y = canopy.xl(2:end); % vertical depth in canopy, 0 - top, -1 - bottom 

%% radiation

figure

subplot(2, 2, 1)

plot(Ps, canopy.xl, 'DisplayName', 'gap.Ps (sunlit)')
hold on
plot(1-Ps, canopy.xl, 'DisplayName', '1 - gap.Ps (shaded)')

legend
ylabel('depth in canopy')
xlabel('probability, -')
title('probability gap fraction', 'per leaf layer')


for p = 2:4
    
    subplot(2, 2, p)

    % aPAR in umol photons m-2 s-1
    name = 'aPAR';
    full_name = {"absorbed PAR", "rad.Pnu, rad.Pnh", "per leaf layer"};
    units = '\mumol photons m-2 s-1';
    sunlit = rad.Pnu;
    shaded = rad.Pnh;
    % fapar canopy.Rntot_PAR / rad.EPAR

    if p == 3
        % net radiation
        name = 'Rn';
        full_name = {"shortwave net radiation", "rad.Rnuc, rad.Rnhc" "per leaf layer"};
        units = 'W m-2';
        sunlit = rad.Rnuc;
        shaded = rad.Rnhc;
    elseif p == 4
        % thermal net radiation
        name = 'Rnt';
        full_name = {"thermal net radiation", "rad.Rnuct, rad.Rnhct", "per leaf layer"};
        units = 'W m-2';
        sunlit = rad.Rnuct;
        shaded = rad.Rnhct;
    end
    
    if options.lite == 0
        sunlit = meanleaf(canopy, sunlit, 'angles');
    end

    plot(sunlit, y, 'o-', 'DisplayName', 'sunlit1d')
    hold on
    plot(shaded, y, 'o-', 'DisplayName', 'shaded')
    
    sunlit_w  = Ps(1:nl) .* sunlit;
    shaded_w = (1-Ps(1:nl)) .* shaded;
    

    % plot(sunlit_w, y, '--', 'DisplayName', 'sunlit1d weighted', 'Color', "#0072BD")
    % hold on
    % plot(shaded_w, y, '--', 'DisplayName', 'shaded weighted', 'Color', 	"#D95319")
    
    plot(sunlit_w + shaded_w, y, 'x-', 'DisplayName', 'sum proflie')
    
    legend
    
    ylabel('depth in canopy')
    xlabel(sprintf('%s, %s', name, units))
    title(full_name)
end

set(findall(gcf,'-property','FontSize'),'FontSize', 14)

%% biochemical output

figure


for p = 1:5
    
    subplot(2, 3, p)

    % Ag (gross leaf photosynthesis) in umol CO2 m-2 s-1
    name = 'Ag';
    full_name = {"gross photosynthesis", "bcu.Ag, bch.Ag", "per leaf layer"};
    units = '\mumol CO_2 m-2 s-1';
    sunlit = bcu.Ag;
    shaded = bch.Ag;

    if p == 2
        name = 'Vcmax';
        full_name = {"maximum carboxylation capactiy", "bcu.Vcmax, bch.Vcmax", "per leaf layer"};
        units = '\mumol CO_2 m-2 s-1';
        sunlit = bcu.Vcmax;
        shaded = bch.Vcmax;
    elseif p == 3
        % eta
        name = 'eta';
        full_name = {"eta = Fs / Fo", "bcu.eta, bch.eta", "per leaf layer"};
        units = '-';
        sunlit = bcu.eta;
        shaded = bch.eta;

    elseif p == 4
        % RH
        name = 'RH_{leaf}';
        full_name = {"relative humidity (RH=0.64)", "bcu.RH, bch.RH", "per leaf layer"};
        units = '-';
        sunlit = bcu.RH;
        shaded = bch.RH;
    elseif p == 5
        % gsw
        name = 'gsw';
        full_name = {"stomatal conductance", "bcu.gs, bch.gs", "per leaf layer"};
        units = 'mol H_2O m-2 s-1';
        sunlit = bcu.gs;
        shaded = bch.gs;
    
    % elseif p == 6
    %     % Ja
    %     name = 'Ja';
    %     full_name = {"ETR actual", "per leaf layer"};
    %     units = '\mumol photons m-2 s-1';
    %     sunlit = bcu.Ja;
    %     shaded = bch.Ja;
    end
    
    if options.lite == 0
        sunlit = meanleaf(canopy, sunlit, 'angles');
    end

    plot(sunlit, y, 'o-', 'DisplayName', 'sunlit1d')
    hold on
    plot(shaded, y, 'o-', 'DisplayName', 'shaded')
    
    sunlit_w  = Ps(1:nl) .* sunlit;
    shaded_w = (1-Ps(1:nl)) .* shaded;
    

    % plot(sunlit_w, y, '--', 'DisplayName', 'sunlit1d weighted', 'Color', "#0072BD")
    % hold on
    % plot(shaded_w, y, '--', 'DisplayName', 'shaded weighted', 'Color', 	"#D95319")
    
    plot(sunlit_w + shaded_w, y, 'x-', 'DisplayName', 'sum proflie')
    
    % legend
    
    ylabel('depth in canopy')
    xlabel(sprintf('%s, %s', name, units))
    title(full_name)
end


%% thermal

subplot(2,3,6)

name = 'T_{leaf}';
full_name = sprintf("temperatures (Ta = %.0f^oC)", meteo.Ta);
units = 'C';
sunlit = thermal.Tcu;
shaded = thermal.Tch;

if options.lite == 0
    sunlit = meanleaf(canopy, sunlit, 'angles');
end


plot(sunlit, y, 'o-', 'DisplayName', 'sunlit1d')
hold on
plot(shaded, y, 'o-', 'DisplayName', 'shaded')


sunlit_w  = Ps(1:nl) .* (sunlit .^ 4);
shaded_w = (1-Ps(1:nl)) .* (shaded .^ 4);

plot((sunlit_w + shaded_w).^ (1/4), y, 'x-', 'DisplayName', 'weighted proflie')

legend

ylabel('depth in canopy')
xlabel(sprintf('%s, %s', name, units))
title({full_name, 'thermal.Tcu, thermal.Tch', 'per leaf layer'})


set(findall(gcf,'-property','FontSize'),'FontSize', 14)

// src/plotting/v7_2_plot_profiles_radiation.m

%% youtube video 7
% https://youtu.be/lhmP78f5tYg?si=O2i3m-WipOSlY8xP

%% plot profiles
% this script works with the structures available in the workspace after
% the SCOPE run

%% radiation
% we just sample top (1) and bottom (end) 
% Emin_ - downward diffuse
% Eplu_ - upwar diffuse

figure

wl = spectral.wlS * 1e-3; % um

subplot(2, 2, 1)

per_wl = rad.Emin_;
% per_wl = rad.Emins_;  % diffuse that used to be direct
% per_wl = rad.Emind_;  % diffuse that came from the Sun

top = per_wl(1, :);
bottom = per_wl(end, :);

plot(wl, top, 'DisplayName', 'top')
hold on
plot(wl, bottom, 'DisplayName', 'bottom')

legend

ylabel('W m-2 \mum-1')
xlabel('wavelength, \mum')
title({'downwelling (E-) radiation', 'rad.Emin_ (from RTMo)', 'diffuse solar and atmosphere'}, ...
    'Interpreter', 'none')

set(gca, 'XScale', 'log')


subplot(2, 2, 2)

per_wl = rad.Eplu_;
% per_wl = rad.Eplus_;  % diffuse that used to be direct
% per_wl = rad.Eplud_;  % diffuse that came from the Sun

top = per_wl(1, :);
bottom = per_wl(end, :);

plot(wl, top, 'DisplayName', 'top')
hold on
plot(wl, bottom, 'DisplayName', 'bottom')

legend

ylabel('W m-2 \mum-1')
xlabel('wavelength, \mum')
title({'upwelling (E+) radiation', 'rad.Eplu_ (from RTMo)', 'reflected solar and atmosphere'}, ...
    'Interpreter', 'none')

set(gca, 'XScale', 'log')



wl = spectral.wlT * 1e-3;

subplot(2, 2, 3)

per_wl = rad.Emint_;

top = per_wl(1, :);
bottom = per_wl(end, :);

plot(wl, top, 'DisplayName', 'top')
hold on
plot(wl, bottom, 'DisplayName', 'bottom')

legend

ylabel('W m-2 \mum-1')
xlabel('wavelength, \mum')
title({'downwelling thermal (E-) radiation', 'rad.Emint_ (from RTMt_sb or RTMt_planck)', 'emitted by canopy'}, ...
    'Interpreter', 'none')

set(gca, 'XScale', 'log')


subplot(2, 2, 4)

per_wl = rad.Eplut_;

top = per_wl(1, :);
bottom = per_wl(end, :);

plot(wl, top, 'DisplayName', 'top')
hold on
plot(wl, bottom, 'DisplayName', 'bottom')

legend

ylabel('W m-2 \mum-1')
xlabel('wavelength, \mum')
title({'upwelling (E+) thermal radiation', 'rad.Eplut_ (from RTMt_sb or RTMt_planck)', 'emitted by canopy'}, ...
    'Interpreter', 'none')

set(gca, 'XScale', 'log')

set(findall(gcf,'-property','FontSize'),'FontSize', 14)

%% per layer

figure

wl = spectral.wlS * 1e-3; % um

for p=1:2
    per_wl = rad.Emin_;
    per_wl(:, spectral.IwlT) = per_wl(:, spectral.IwlT) + rad.Emint_;
    what = {'downwelling (E-) radiation', 'rad.Emin_ + rad.Emint_'};
    
    if p == 2
        per_wl = rad.Eplu_;
        per_wl(:, spectral.IwlT) = per_wl(:, spectral.IwlT) + rad.Eplut_;
        what = {'upwelling (E+) radiation', 'rad.Eplu_ + rad.Eplut_'};
    end
    

    subplot(1, 2, p)
    plot(wl, per_wl(1, :), 'DisplayName', 'top')
    hold on
    for i=1:(nl - 1)
        if rem(i, 5) == 0  % each 5th layer
            plot(wl, per_wl(i, :), 'DisplayName', sprintf('layer %d', i))
        end
    
    end
    plot(wl, per_wl(end, :), 'DisplayName', 'bottom')
    
    legend
    set(gca, 'XScale', 'log')
    title(what, 'Interpreter', 'none')
    
    ylabel('W m-2 \mum-1')
    xlabel('wavelength, \mum')
end


set(findall(gcf,'-property','FontSize'),'FontSize', 14)


%% components

figure

wl = spectral.wlS * 1e-3; % um

subplot(2, 2, 1)

per_wl = rad.Emin_;
% per_wl = rad.Emins_;  % diffuse that used to be direct
% per_wl = rad.Emind_;  % diffuse that came from the Sun

top = per_wl(1, :);
bottom = per_wl(end, :);

plot(wl, top, 'DisplayName', 'top')
hold on
plot(wl, bottom, 'DisplayName', 'bottom')

legend

ylabel('W m-2 \mum-1')
xlabel('wavelength, \mum')
title({'downwelling (E-) radiation', 'rad.Emin_ (from RTMo)', 'diffuse solar and atmosphere'}, ...
    'Interpreter', 'none')

set(gca, 'XScale', 'log')


subplot(2, 2, 2)

per_wl = rad.Eplu_;
% per_wl = rad.Eplus_;  % diffuse that used to be direct
% per_wl = rad.Eplud_;  % diffuse that came from the Sun

top = per_wl(1, :);
bottom = per_wl(end, :);

plot(wl, top, 'DisplayName', 'top')
hold on
plot(wl, bottom, 'DisplayName', 'bottom')

legend

ylabel('W m-2 \mum-1')
xlabel('wavelength, \mum')
title({'upwelling (E+) radiation', 'rad.Eplu_ (from RTMo)', 'reflected solar and atmosphere'}, ...
    'Interpreter', 'none')

set(gca, 'XScale', 'log')

wl = spectral.wlS * 1e-3; % um


subplot(2, 2, 3)

per_wl = rad.Emin_;
% per_wl = rad.Emins_;  % diffuse that used to be direct
% per_wl = rad.Emind_;  % diffuse that came from the Sun

top = per_wl(1, :);
bottom = per_wl(end, :);



% plot(wl, top)
hold on
plot(wl, rad.Esun_, 'Color', "#EDB120")
plot(wl, rad.Esky_, 'Color', "#0072BD")

stacked = [rad.Emins_(end,:); rad.Emind_(end,:)]; 
a = area(wl, stacked');
a(1).FaceColor = "#EDB120";
a(2).FaceColor = "#0072BD";

set(gca, 'XScale', 'log')

legend({'rad.Esun_ (direct solar light)', 'rad.Esky_ (diffuse solar light)',...
    ['rad.Eplus_' newline 'diffuse that used to be direct'], ...
    ['rad.Emind_' newline 'diffuse that came from the Sun']}, 'Interpreter', 'none')

title({'components of downwelling (E-) radiation [bottom]', 'rad.Emin_ (from RTMo)', 'diffuse solar and atmosphere'}, ...
    'Interpreter', 'none')
ylabel('W m-2 \mum-1')

subplot(2, 2, 4)

per_wl = rad.Eplu_;
% per_wl = rad.Eplus_;  % diffuse that used to be direct
% per_wl = rad.Eplud_;  % diffuse that came from the Sun

top = per_wl(1, :);
bottom = per_wl(end, :);

stacked = [rad.Eplus_(1,:); rad.Eplud_(1,:)]; 

% plot(wl, top)
% hold on
a = area(wl, stacked');
a(1).FaceColor = "#EDB120";
a(2).FaceColor = "#0072BD";

set(gca, 'XScale', 'log')

legend({['rad.Eplus_' newline 'diffuse that used to be direct'], ...
    ['rad.Eplud_' newline 'diffuse that came from the Sun']}, 'Interpreter', 'none')

title({'components of upwelling (E+) radiation [top]', 'rad.Eplu_ (from RTMo)', 'reflected solar and atmosphere'}, ...
    'Interpreter', 'none')
ylabel('W m-2 \mum-1')

set(findall(gcf,'-property','FontSize'),'FontSize', 14)

// src/plotting/v8_plot_mSCOPE.m

%% youtube video 8
% https://youtu.be/f543CkCSIUg?si=etRXbAaf7dUtxtX1

%% mSCOPE compute different reflectance for leaf layers

figure

wl = spectral.wlS * 1e-3; 

for i=1:canopy.nlayers
    plot(wl, leafopt.refl(i, :), 'DisplayName', 'leaf reflectance')
    hold on
end


set(gca, 'XScale', 'log')


xlabel('wavelength, \mum')
ylabel('reflectance, -')
title('leaf reflectance', 'per layer')

set(findall(gcf,'-property','FontSize'),'FontSize', 14)

%%

nl = canopy.nlayers;
Ps = gap.Ps;

y = canopy.xl(2:end); % vertical depth in canopy, 0 - top, -1 - bottom 

%% radiation

figure

subplot(2, 2, 1)

plot(Ps, canopy.xl, 'DisplayName', 'gap.Ps (sunlit)')
hold on
plot(1-Ps, canopy.xl, 'DisplayName', '1 - gap.Ps (shaded)')

legend
ylabel('depth in canopy')
xlabel('probability, -')
title('probability gap fraction', 'per leaf layer')



for p = 2:4
    
    subplot(2, 2, p)

    % aPAR in umol photons m-2 s-1
    name = 'aPAR';
    full_name = {"absorbed PAR", "rad.Pnu, rad.Pnh", "per leaf layer"};
    units = '\mumol photons m-2 s-1';
    sunlit = rad.Pnu;
    shaded = rad.Pnh;
    % fapar canopy.Rntot_PAR / rad.EPAR


    if p == 3
        name = 'Ag';
        full_name =  {"gross photosynthesis", "bcu.Ag, bch.Ag", "per leaf layer"};
        units = '\mumol CO_2 m-2 s-1';
        sunlit = bcu.Ag;
        shaded = bch.Ag;
    elseif p == 4
        % eta
        name = 'eta';
        full_name = {"eta = Fs / Fo", "bcu.eta, bch.eta", "per leaf layer"};
        units = '-';
        sunlit = bcu.eta;
        shaded = bch.eta;
    end
    
    if options.lite == 0
        sunlit = meanleaf(canopy, sunlit, 'angles');
    end

    plot(sunlit, y, 'o-', 'DisplayName', 'sunlit1d')
    hold on
    plot(shaded, y, 'o-', 'DisplayName', 'shaded')
    
    sunlit_w  = Ps(1:nl) .* sunlit;
    shaded_w = (1-Ps(1:nl)) .* shaded;
    

    % plot(sunlit_w, y, '--', 'DisplayName', 'sunlit1d weighted', 'Color', "#0072BD")
    % hold on
    % plot(shaded_w, y, '--', 'DisplayName', 'shaded weighted', 'Color', 	"#D95319")
    
    plot(sunlit_w + shaded_w, y, 'x-', 'DisplayName', 'sum proflie')
    
    legend
    
    ylabel('depth in canopy')
    xlabel(sprintf('%s, %s', name, units))
    title(full_name)
end

set(findall(gcf,'-property','FontSize'),'FontSize', 14)

// src/plotting/v9_plot_sif.m

%% youtube video 9
% https://youtu.be/CpE07VEjs4Q?si=r46ltIpa94ljWff_

%% fluorescence background
% emitted 
% (1) at photosystem level (inside chloroplasts)
% (2) in all directions (hemispherical)
% reabsorbed 
% (1) within the leaf it was emitted
% (2) within the canopy (by other leaves)
% (3) by soil


%% hemispherical (theoretical + real)

figure
wl = spectral.wlF;  % nm

% patch([x fliplr(x)], [y1 fliplr(y2)], 'g')
area(wl, rad.EoutFrc_, 'DisplayName', ['rad.EoutFrc_', newline ...
    'photosystem level (without within-leaf reabsorption)' newline ...
    'fluorescence_ReabsCorr.csv'])
hold on

area(wl, rad.Femleaves_, 'DisplayName', ['rad.Femleaves_', newline ...
    'all leaves forward+backward (without canopy/soil reabsorption)' newline ...
    'fluorescence_AllLeaves.csv'])

area(wl, rad.EoutF_, 'DisplayName', ['rad.EoutF_', newline ...
    'SIF hemispherical (at the top of canopy)' newline...
    'fluorescence_hemis.csv'])


legend('Interpreter','none')


xlabel('wavelength, nm')
ylabel('E, mW m-2 nm-1')  % == W m-2 um-1

set(findall(gcf,'-property','FontSize'),'FontSize', 14)

title({'SIF', 'hemispherically integrated'})

%% stacked: in observation direction

figure
wl = spectral.wlF;  % nm

plot(wl, rad.LoF_, 'o', 'DisplayName', 'in observation direction')
hold on

stacked = [rad.LoF_sunlit, rad.LoF_shaded, rad.LoF_scattered, rad.LoF_soil];
area(wl, stacked)

legend({['rad.LoF_' newline 'fluorescence.csv'],...
    'rad.LoF_sunlit', 'rad.LoF_shaded', 'rad.LoF_scattered', 'rad.LoF_soil'})

legend('Interpreter','none')

xlabel('wavelength, nm')
ylabel('L, mW m-2 nm-1 sr-1')  % == W m-2 um-1

set(findall(gcf,'-property','FontSize'),'FontSize', 14)

title({'SIF', 'in the observation direction'})


%% escape probability sigmaF

figure

subplot(2, 2, 1)

wl = spectral.wlF;  % nm

plot(wl, rad.LoF_ * pi, 'o')


xlabel('wavelength, nm')
ylabel('L * pi, mW m-2 nm-1')

title({'SIF in the observation direction (measured)', 'rad.LoF_ multiplied by pi'}, 'Interpreter', 'none')

%
subplot(2, 2, 3)

plot(wl, rad.sigmaF)

xlabel('wavelength, nm')
ylabel('sigmaF, -')
title({'SIF escape probability (rad.sigmaF)', 'escape from photosystem to sensor'}, 'Interpreter', 'none')

% plot(wl, rad.EoutFrc_ .* rad.sigmaF' / pi)
% hold on
% plot(wl,  rad.LoF_)


%%%%% reabsorption corrected

subplot(2, 2, [2, 4])


% from ETR and fqe of PSII
plot(wl, rad.EoutFrc_, 'x')


xlabel('wavelength, nm')
ylabel('E, mW m-2 nm-1')
title({'SIF hemispherically integrated (theoretical)', ...
    'rad.EoutFrc_ = rad.LoF_ * pi / rad.sigmaF'}, 'Interpreter', 'none')


set(findall(gcf,'-property','FontSize'),'FontSize', 14)

%% apparent reflectance and radiance

figure

subplot(2,2,1)
wl = spectral.wlS * 1e-3; 

plot(wl, rad.refl, 'DisplayName', 'rad.refl (reflectance)')
hold on
plot(wl, rad.reflapp, 'DisplayName', 'rad.reflapp (apparent reflectance)')

set(gca, 'XScale', 'log')

title('canopy reflectances')
xlabel('wavelength, \mum')
ylabel('reflectance, -')
legend('Interpreter', 'none')


subplot(2,2,2)
plot(wl, rad.reflapp - rad.refl, 'Color', "#EDB120")

set(gca, 'XScale', 'log')

title('apparent minus true reflectance', 'rad.reflapp - rad.refl', 'Interpreter','none')
xlabel('wavelength, \mum')
ylabel('reflectance, -')


subplot(2,2,3)

plot(wl, rad.Lotot_, 'DisplayName', 'rad.Lotot_ (radiance)')
hold on
plot(wl, rad.Lototf_, 'DisplayName', 'rad.Lototf_ (apparent radiance)')

set(gca, 'XScale', 'log')

title('canopy radiance in observation direction')
xlabel('wavelength, \mum')
ylabel('L, W m-2 \mum-1 sr-1')
legend('Interpreter', 'none')

subplot(2,2,4)

plot(wl, rad.Lototf_ - rad.Lotot_, 'Color', "#EDB120")

set(gca, 'XScale', 'log')

title('apparent minus true radiance', 'rad.Lototf_ - rad.Lotot_', 'Interpreter','none')
xlabel('wavelength, \mum')
ylabel('L, W m-2 \mum-1 sr-1')


set(findall(gcf,'-property','FontSize'),'FontSize', 14)


for i=1:4
    subplot(2,2,i)
    xlim([640, 850]*1e-3)
end

// src/supporting/Monin_Obukhov.m

function L = Monin_Obukhov(constants,meteo,H)

L           = -constants.rhoa*constants.cp*meteo.ustar.^3.* ...
            (meteo.Ta+273.15)./(constants.kappa*constants.g*H);           % [1]
     
L(isnan(L)) = -1E6;

// src/supporting/Planck.m

function Lb = Planck(wl,Tb,em)

    c1 = 1.191066e-22;
    c2 = 14388.33;
    if nargin<3
        em = ones(size(Tb));
    end
    
    Lb = em.* c1*(wl*1e-9).^(-5)./(exp(c2./(wl*1e-3*Tb))-1);
    
end    

// src/supporting/Sint.m


function int = Sint(y,x)

    % Simpson integration
    % x and y must be any vectors (rows, columns), but of the same length
    % x must be a monotonically increasing series
    
    % WV Jan. 2013, for SCOPE 1.40
    
    nx   = length(x);
    if size(x,1) == 1
        x = x';
    end
    if size(y,1) == size(x,1)
        y = y';
    end

    step = x(2:nx) - x(1:nx-1);
    mean = .5 * (y(:,1:nx-1) + y(:,2:nx));
    int  = mean * step;
end

// src/supporting/Soil_Inertia0.m

function [GAM]  = Soil_Inertia0(cs,rhos,lambdas)
% soil thermal inertia 
GAM             = sqrt(cs*rhos*lambdas);                            % soil thermal intertia

// src/supporting/Soil_Inertia1.m

function [GAM] = Soil_Inertia1(SMC)

%soil inertia method by Murray and Verhoef (

%% parameters

theta_s = 0.435; %(saturated water content, m3/m3)
Sr = SMC/theta_s;

%fss = 0.58; %(sand fraction)
gamma_s = 0.96; %(soil texture dependent parameter)
dels = 1.33; %(shape parameter)


ke = exp(gamma_s*(1- power(Sr,(gamma_s - dels))));

phis  = 0.435; %(phis == theta_s)
lambda_d = -0.56*phis + 0.51;

QC = 0.60; %(quartz content)
lambda_qc = 7.7;  %(thermal conductivity of quartz, constant)

lambda_s = (lambda_qc^(QC))*lambda_d^(1-QC);
lambda_wtr = 0.57;   %(thermal conductivity of water, W/m.K, constant)

lambda_w = (lambda_s^(1-phis))*lambda_wtr^(phis);

lambdas = ke*(lambda_w - lambda_d) + lambda_d;

Hcs = 2*10^6;
Hcw = 4.2*10^6;

Hc = (Hcw * SMC)+ (1-theta_s)*Hcs;

GAM = sqrt(lambdas.*Hc);

// src/supporting/aggreg.m

function [M] = aggreg(atmfile,SCOPEspec)

% Aggregate MODTRAN data over SCOPE bands by averaging (over rectangular band
% passes)

% Read .atm file with MODTRAN data
s   = importdata(atmfile);
wlM = s.data(:,2);
T   = s.data(:,3:20);

% Extract 6 relevant columns from T

%  1: <Eso*costts/pi>
%  3: <rdd>
%  4: <tss>
%  5: <tsd>
% 12: <tssrdd>
% 16: <La(b)>

U     = [T(:,1) T(:,3) T(:,4) T(:,5) T(:,12) T(:,16)];

nwM   = length(wlM);

nreg  = SCOPEspec.nreg;
streg = SCOPEspec.start;
enreg = SCOPEspec.end;
width = SCOPEspec.res;

% Nr. of bands in each region

nwreg = int32((enreg-streg)./width)+1;

off   = int32(zeros(nreg,1));

for i=2:nreg
    off(i) = off(i-1)+nwreg(i-1);
end

nwS = sum(nwreg);
n   = zeros(nwS,1);    % Count of MODTRAN data contributing to a band
S   = zeros(nwS,6);    % Intialize sums

%k   = int32(0);
j   = int32(zeros(nreg,1));  % Band index within regions

for iwl = 1:nwM
    w   = wlM(iwl);    % MODTRAN wavelength
    for r = 1:nreg
        j(r) = int32(round(w-streg(r))./(width(r)))+1;
        if j(r)>0 && j(r)<=nwreg(r)                 % test if index is in valid range
            k      = j(r)+off(r);                   % SCOPE band index
            S(k,:) = S(k,:)+U(iwl,:);               % Accumulate from contributing MODTRAN data
            n(k)   = n(k)+1;                        % Increment count
        end
    end
end

M = zeros(size(S,1),6);
for i = 1:6
    M(:,i) = S(:,i)./n;      % Calculate averages per SCOPE band
end

end

// src/supporting/calc_brdf.m

function directional = calc_brdf(constants,options,directional,spectral,angles,atmo,soil,leafopt,canopy,meteo,thermal,bcu,bch)

% simulates observations from a large number of viewing angles
% modified: 30 April 2020, CvdT, removed repeated angle combinations.

%% input
tts             = angles.tts;
psi_hot         = [0    ; 0     ;0     ;0     ;0    ;2  ;358];  % [noa_o]           angles for hotspot oversampling
tto_hot         = [tts  ; tts+02;tts+04;tts-02;tts-4;tts;tts];  % [noa_o]           angles for hotspot oversampling

psi_plane       = [000*ones(6,1);180*ones(6,1);090*ones(6,1);270*ones(6,1)];%       angles for plane oversampling
tto_plane       = [10:10:60     , 10:10:60    , 10:10:60    , 10:10:60]';   %       angles for plane oversampling

psi             = [directional.psi; psi_hot; psi_plane];
tto             = [directional.tto; tto_hot; tto_plane];

[~,u]           = unique([psi tto],'rows');
directional.psi = psi(u);
directional.tto = tto(u);
na              = length(u);

%% allocate memory
directional.brdf_       = zeros(length(spectral.wlS),na);      % [nwlS, no of angles]  
directional.Eoutte      = zeros(1,na);                         % [1, no of angles]  
directional.BrightnessT = zeros(1,na);                         % [1, no of angles] 
directional.LoF_        = zeros(length(spectral.wlF),na);      % [nwlF, no of angles] 
directional.Lot_        = zeros(length(spectral.wlT),na);      % [nwlF, no of angles] 

%% other preparations
directional_angles = angles;

%% loop over the angles
for j=1:na
    
    %optical BRDF
    directional_angles.tto = directional.tto(j);
    directional_angles.psi = directional.psi(j);
    [directional_rad,directional_gap] = RTMo(spectral,atmo,soil,leafopt,canopy,directional_angles,constants,meteo,options);
    directional.refl_(:,j)  = directional_rad.refl;         % [nwl]            reflectance (spectral) (nm-1)
    directional.rso_(:,j)  = directional_rad.rso;           % [nwl]            BRDF (spectral) (nm-1)
    
    % thermal directional brightness temperatures (Planck)
    
    if options.calc_planck
        directional_rad = RTMt_planck(spectral,directional_rad,...
            soil,leafopt,canopy,directional_gap,...
            thermal.Tcu,thermal.Tch,thermal.Tsu,thermal.Tsh);
        directional.Lot_(:,j)      = directional_rad.Lot_(spectral.IwlT) + directional_rad.Lo_(spectral.IwlT);        % [nwlt]           emitted plus reflected diffuse     radiance at top
    end
    
    if options.calc_fluor      
        directional_rad = RTMf(constants,spectral,directional_rad,soil,leafopt,canopy,directional_gap,directional_angles,bcu.eta,bch.eta);    
        directional.LoF_(:,j)              = directional_rad.LoF_;
    end % {if calcfluor}
    if options.calc_xanthophyllabs
        directional_rad = RTMz(constants,spectral,directional_rad,soil,leafopt,canopy,directional_gap,directional_angles,bcu.Kn,bch.Kn);
        directional.Lo_ = directional_rad.Lo_;
    end
                
end % {for angles}

// src/supporting/calc_rssrbs.m

function [rss,rbs] = calc_rssrbs(SMC,LAI,rbs)

rss            = 11.2*exp(42*(0.22-SMC));
rbs            = rbs*LAI/3.3;

// src/supporting/calczenithangle.m

function [Fi_s,Fi_gs,Fi_g,Omega_s] = calczenithangle(Doy,t,Omega_g,Fi_gm,Long,Lat)
%
% author: Christiaan van der Tol (c.vandertol@utwente.nl)
% date:     Jan 2003
% update:   Oct 2008 by Joris Timmermans (j_timmermans@itc.nl): 
%               - corrected equation of time
%           Oct 2012 (CvdT) comment: input time is GMT, not local time!
%
% function [Fi_s,Fi_gs,Fi_g]= calczenithangle(Doy,t,Omega_g,Fi_gm,Long,Lat)
%
% calculates pi/2-the angle of the sun with the slope of the surface.
%
% input:
% Doy       day of the year
% t         time of the day (hours, GMT)
% Omega_g   slope azimuth angle (deg)
% Fi_gm     slope of the surface (deg)
% Long      Longitude (decimal)
% Lat       Latitude (decimal)
%
% output:
% Fi_s      'classic' zenith angle: perpendicular to horizontal plane
% Fi_gs     solar angle perpendicular to surface slope
% Fi_g      projected slope of the surface in the plane through the solar beam and the vertical
% 

%parameters (if not already supplied)
if nargin<6
    Long        =   13.75;                      % longitude
    Lat         =   45.5;                       % latitude
    if (nargin<4)
        Omega_g =   210;                        % aspect
        Fi_gm   =   30;                         % slope angle
    end
end

%convert angles into radials
G               =   (Doy-1)/365*2*pi;           % converts day of year to radials
Omega_g         =   Omega_g/180*pi;             % converts direction of slope to radials
Fi_gm           =   Fi_gm/180*pi;               % converts maximum slope to radials
Lat             =   Lat/180*pi;                 % converts latitude to radials

%computes the declination of the sun
d               =   0.006918-0.399912*cos(G  )+ 0.070247*sin(G  )- ...
                     0.006758*cos(2*G)+ 0.000907*sin(2*G)- ...
                     0.002697*cos(3*G)+ 0.00148*sin(3*G);
                                
%Equation of time
Et              =   0.017 + .4281 * cos(G) - 7.351 * sin(G) - 3.349 * cos(2*G) - 9.731 * sin(2*G);

%computes the time of the day when the sun reaches its highest angle                                
tm              =   12+(4*(-Long)-Et)/60;      % de Pury and Farquhar (1997), Iqbal (1983)

%computes the hour angle of the sun
Omega_s         =   (t-tm)/12*pi;

%computes the zenithangle (equation 3.28 in De Bruin)
Fi_s            =   acos(sin(d)*sin(Lat)+cos(d)*cos(Lat).*cos(Omega_s));

%computes the slope of the surface Fi_g in the same plane as the solar beam
Fi_g            =   atan(tan(Fi_gm).*cos(Omega_s-Omega_g));

%computes the angle of the sun with the vector perpendicular to the surface
Fi_gs           =   Fi_s + Fi_g;

// src/supporting/count_k.m

function [vnew]=count_k(nvars,v,vmax,id)

% nvars = number of digits
% v     = current vector of digits
% vmax  = maximum values of digits
% id    = starting digit
% vnew  = new vector of digits

i=id;

% starting at id, set digits which are at its maximum equal to 1
% first digit that is not at its maximum is incremented

while v(i)==vmax(i)
    v(i)=1;
    i=rem(i,nvars)+1;
end
v(i)=rem(v(i),vmax(i))+1;
vnew=v;
end

// src/supporting/e2phot.m

function molphotons = e2phot(lambda,E,constants)
%molphotons = e2phot(lambda,E) calculates the number of moles of photons
%corresponding to E Joules of energy of wavelength lambda (m)

e           = ephoton(lambda,constants);
photons     = E./e;
molphotons  = photons./constants.A;
return;

// src/supporting/ephoton.m

function E = ephoton(lambda,constants)
%E = phot2e(lambda) calculates the energy content (J) of 1 photon of
%wavelength lambda (m)

h       = constants.h;           % [J s]         Planck's constant
c       = constants.c;           % [m s-1]       speed of light
E       = h*c./lambda;           % [J]           energy of 1 photon
return;

// src/supporting/fixedp_brent_ari.m

function [b, err2, fcounter] = fixedp_brent_ari(func, x0, corner, tolFn, verbose)
% Find a fixed point of func(x) using Brent's method, as described by Brent 1971
%   func is a single-argument function, f(x) that returns a value the same size as x:
%        The goal is to find f(x) = x (or for Brent, f(x) - x = 0).
%   x0 is the initial guess (or 2 x n matrix if we want to generalize)
%   tol is the tolerance in x (or if two-valued, x, f(x)? )
%   corner (optional) is a known "edge" in the function that could slow down the algorithm
%      if specified and the first two points include the corner, the corner will be substituted as a starting point.
%
%  Written by: Ari Kornfeld, 2016-10

tolx_1 = 0; % we don't want to converge in x
subsetLimit = 0; % never use subsets of x
accelBisection = false; % should we divide the interval by ever-increasing fractions (1/2, 1/3, 1/4) when consecutive calls make no improvement?
accelLimit = 3; % but reset if more than accelLimit bisections don't help.
rotatePrev = true; % do we use Brent's original method of wiping out "xprev" in certain instances (forcing secant) or alway preserve a third value (allowing more inverse-quad).
iter_limit = 100;

if nargin < 5
    verbose = false;
end

if nargin < 4
    tolFn = eps; % default to FP precision: 2^-52 = 2.2204e-16
end
if nargin < 3
    corner = []; % default, no corners
end
track_fcount = (nargout > 2) || verbose;

recompute_b = false; % ensure that we return after a call on func(s); this allows us to skip the re-call in the main body (since func sets globals)

%  We keep track of the best value of x (b), the past two iterations (c, d),
%    and one contrapoint (a), i.e. a point that stradles the fixed point (i.e. err1 has the opposite sign of err2)
%  err1, err2, ... are the corresponding "y" values for y = f(x) - x (i.e. distance from the fixed point)
% We start by finding bracketing points a, b
a = x0;
[err1, b] = func(a);  % guess the second point by computing x_1 = f(x_0)

% special case: func may return a vector even though 'a' was scalar (because other variables not visible here were nonscalar).
% If so, expand a:
if length(a) == 1 && length(b) > 1  % anything else is an error
    a = repmat(a, size(b));  % or a * ones(size(b))?
end

err2 = func(b);
err2(isnan(err2)) = 0;  % if isnan, we're done, i.e. we shouldn't try
if track_fcount
    % count each item separately: count only the "necessary" iterations
    fcounter = 2*ones(size(b));
end

err_outside_tol = abs(err2) > tolFn;
if ~any(err_outside_tol)
    return  % we're done!
end
% ELSE
recompute_b = true;  % we'll be messing with it

% Now confirm that the two first guesses bracket zero.
%  NOTE: the algorithm may still succeed w/o bracketting, though it's not guaranteed.
not_bracketting_zero = (sign(err1) == sign(err2)) & err_outside_tol;
if any(not_bracketting_zero)
    %warning( 'Not all initial guesses bracket zero. Will fix it now.' );
    
    % first try a simple secant extrapolation
    x1 = b - err2.*(b - a)./(err2 - err1);
    err_x1 = func(x1);
    if track_fcount
        fcounter = fcounter + not_bracketting_zero; % count only the ones that needed this fcall
    end
    
    % since sign(err1) == sign(err2), compare the new value to either of those
    use_x1 = (sign(err_x1) ~= sign(err1)) & not_bracketting_zero;
    if any(use_x1)
        % save the better of the two original points into 'a'
        swap_to_a = (abs(err2) < abs(err1)  & use_x1);
        a(swap_to_a) = b(swap_to_a);      err1(swap_to_a) = err2(swap_to_a);
        % then put the new contrapoint into 'b'
        b(use_x1) = x1(use_x1);           err2(use_x1) = err_x1(use_x1);
    end
    
    % recompute a_too_high and iterate if necessary
    err_outside_tol = min(abs(err1), abs(err2)) > tolFn;
    not_bracketting_zero = (sign(err1) == sign(err2)) & err_outside_tol;
   % make 'a' the lower value, to make the rest simpler
   if any(not_bracketting_zero)
       swap_to_a = (err2 < err1  & not_bracketting_zero);
       [a(swap_to_a), b(swap_to_a)] = deal(b(swap_to_a), a(swap_to_a));
       [err1(swap_to_a), err2(swap_to_a)] = deal( err2(swap_to_a), err1(swap_to_a) );
   end
   % if both values > 0, need to find a negative value:
   both_positive = err1 > 0 & not_bracketting_zero;
   ntries=1; % up to 10 tries
   while any(both_positive)
        % err1 < err2, so assuming a monotonic fn, we can decrease err1 by increasing distance from 'b'
       diffab = b(both_positive) - a(both_positive);  
       a(both_positive) = a(both_positive) - diffab; % walk out in steps that double with each iteration
       % recompute the new a values (note, it might be smarter to shift a[n-1] into b when we're done
       err1 = func(a);
       if track_fcount
           fcounter = fcounter + both_positive;
       end
       
       % for a severely ill-behaved function err1 can go in the wrong direction as we move apart, so fix it now
       swap_to_a = (err2 < err1  & not_bracketting_zero);
       [a(swap_to_a), b(swap_to_a)] = deal(b(swap_to_a), a(swap_to_a));
       [err1(swap_to_a), err2(swap_to_a)] = deal( err2(swap_to_a), err1(swap_to_a) );
       err_outside_tol = min(abs(err1), abs(err2)) > tolFn;
       not_bracketting_zero = (sign(err1) == sign(err2) & err_outside_tol);
       both_positive = not_bracketting_zero;
       if any(both_positive) && ntries > 10
            error('Couldn''t find contrapoint in 10 tries!')
        end
        ntries = ntries + 1;
   end
    
%    ntries=1; % if using while loop for b
    both_negative = err2 < 0 & not_bracketting_zero;
    if any(both_negative)
        %  for Ci, A(0) -> Ci >= 0 so  f(0) - 0 >= 0
        b(both_negative) = 0; % just go to zero and be done NOT GENERAL!!!!
        err2 = func(b);  % for now don't calculate on subsets
        if track_fcount
            fcounter = fcounter + both_negative;
        end
    end
    
    recompute_b = true; % we can no longer be certain that s (b) is the best)
end

if ~isempty(corner)
    % special case: guesses that bracket a "corner" may result in very slow convergence
    %  if we replace an endpoint with corner, the remaining iterations should be much simplified
    bracket_corner = sign(corner - a) ~= sign(corner - b); %( a < corner & b > corner) | ( a > corner & b < corner);
    if any(bracket_corner)
        % replace one endpoint with corner:
        x1 = b; % initialize it for the function call; we'll only use the bracket_corner elements
        if length(corner) == 1 % as above, make sure value is same size as x
            corner = repmat(corner, size(b));
        end
        x1(bracket_corner) = corner(bracket_corner);
        
        errCorner = func(x1);
        if track_fcount
            fcounter = fcounter + bracket_corner;
        end
        
        % sort the result into a b so that the sign is preserved
        save_into_b = bracket_corner & (sign(errCorner) == sign(err2)); % error on corner and b have same sign
        save_into_a = bracket_corner & ~save_into_b;
        [a(save_into_a),b(save_into_b)] = deal(x1(save_into_a), x1(save_into_b));
        [err1(save_into_a),err2(save_into_b)] = deal(errCorner(save_into_a), errCorner(save_into_b));

        recompute_b = true; % we can no longer be certain that s (i.e. matches the last call to computeA)
    end
end


% Brent tolerance:
% BRENT: 2 * eps * abs(b)+ toler %NUMERICAL RECIPES: 2 * eps * abs(b)+ 0.5*toler
tolx = 2* max(1, abs(b)).*tolx_1; % MATLAB's criterea (and also Wilkins 2013, except for the factor of 4)
err_outside_tol =  0.5 .* abs(a - b) > tolx & min(abs(err1), abs(err2)) > tolFn ;


% make sure that 'b' is the best 
err1_is_best = abs(err2) > abs(err1); % and therefore needs to be swapped
if any(err1_is_best)
    [a(err1_is_best), b(err1_is_best)] = deal(b(err1_is_best), a(err1_is_best));
    [err1(err1_is_best), err2(err1_is_best)] = deal(err2(err1_is_best), err1(err1_is_best));
    recompute_b = true;
end

% initialize the search array with the current best guess
%[s, err_s] = deal(b, err2); 
ab_gap =  (a - b);  % do this only after we've sorted a and b

% since b is the current best guess, 'a' therefore stands in as the "previous best"
[c, err3] = deal(a, err1);  % the "previous" value of b
best_is_unchanged = abs(err2) == abs(err1);

% initialize additional vectors to nan
%[p, q] = deal(nan(size(b)));  % components of the new increment
[xstep, xstep1] = deal(3*ab_gap); % the step_size of one- or two-steps ago; initialize to prevent triggering on the first round
q = ones(size(b)); % needed for first iteration; don't use nan: we need p/q to be a valid number
p = 0 .* q;
%--------------------------------------------------------
%  MAIN LOOP
%  note: Stage 2013 says we should add   abs(a - b)<= 0.5(a+b)*eps as a stopping condition
%   to avoid machine accuracy issues when abs(x_root) >> 1.
%     fzero.m uses: 0.5*abs(a-b) <= 2.0*tol*max(abs(b),1.0) [adjusting for the fact that c in fzero is our a]
%used_bisection = true(size(b));  % 'mflag': did we use bisection on the previous round (starts true)
counter = 0;
accel_bi = zeros(size(b));
if verbose
    fprintf('init''l ')
    fprintf('%d: a: %9.5g (%9.3g), b:  %9.5g (%9.3g),c: %9.5g (%9.3g), s: %9.3g\n', fcounter, a, err1, b, err2, c, err3, err2); % a, b, c, s);
end
while any( err_outside_tol )
    % 0. Setup
    %*** NOTE: See 2013 Wilkins about an objection to the magnitude of tol in these tests 
    xstep2 = xstep1; % Advance the record of previous step-sizes; Brent calls this "e"
    xstep1 = xstep;
    
    %ab_gap =  (a - b); % done at end of loop since it's needed for the exit test
    p = 0.*p; % clear p, xstep (this is a bit faster than zeros() or even p(:)=0 in R2015b )
    xstep = 0.*xstep;
%    [p, xstep] = deal(zeros(size(b))); % don't use nan, so p/q is a valid number
    %q = ones(size(b)); % don't worry about q: if not used, p/q = 0
    
    %use_bisection = (abs(xstep2) < tol) | (abs_err_s >= abs_err_s1); % if the latest err is larger than the previous iteration give up on interpolation
    %abs_err_s1 = abs_err_s;
    use_bisection = (abs(xstep2) < tolx) | best_is_unchanged ; % if the latest err is larger than the previous iteration give up on interpolation
    
    r2 = err2 ./ err1; % f(b)/f(c) in Brent - see my comments in secant
    
    % the new guess will be stored in 's' using either secant, inverse quadratic, or bisection:
    try_interp = ~use_bisection & err_outside_tol;
    %err3_close_enough = abs(ab_gap)*20 >= abs(err3 - err1); % if 3rd best is too far, just use two 
    quad_is_safe = (err1 ~= err3 & err2 ~= err3); % (note err1 never equals err2 -- since they have opp. signs -- unless they're zero, in which case try_interp is false)

    %-----
    % 1. Inverse quadratic Method
    % when three values are available, use inverse quadratic method (note: fzero has a more clever algorithm?)
    use_quad = try_interp & quad_is_safe; % see prev note about (err1 ~= err2) 
    if any(use_quad)
        % this way is 3x faster than with subsetting!
           % (defining intermediate difference variables doesn't make much if any speed difference.)
%         s1 =  a .* err2 .* err3 ./ ((err1 - err2).*(err1 - err3)) ...
%              + b .* err1 .* err3 ./ ((err2 - err1).*(err2 - err3))  ...
%              + c .* err1 .* err2 ./ ((err3 - err1).*(err3 - err2));
%         s(use_quad) = s1(use_quad);
        r1 = err3 ./ err1; % Brent Section 4, not ALGOL, notation, swapping a/c from Brent's notation
        r3 = err2 ./ err3;
        p = r3 .* ( ab_gap .* r1 .* (r1 - r2) - (b - c).*(r2 - 1) );
        q = (r1 - 1) .* (r2 - 1) .* (r3 - 1);
%         p(~use_quad) = 0; % so xstep = 0
%         q(~use_quad) = 1; % so xstep is not nan
%         p(use_quad) = p1(use_quad); % don't bother; we'll overwrite what needs overwriting
%         q(use_quad) = q1(use_quad);
        if verbose, fprintf('**quad '), end
    end
    
    %-----
    % 2. Secant method
    % I've found no case in which doing secant helps when inv. quad failed (it actually makes things worse)
    %s_test = (quad_is_safe & abs(pq1) >= abs(ab_gap) ) ; % if inverse-quad went too far
    use_secant = try_interp & ~quad_is_safe;
    % secant:  b - f(b) * (b - a) / ( f(b) - (fa) ): derivation: using point-slope of a line solve for y=0:
    %      (y - y1) = m(x - x1):  m = (y2-y1)/(x2-x1);  b = y1 - m x1 = y2 - m x2;  
    %      y = 0 => -y1 = m x - m x1 => -mx = y1 - m x1 => x = -y1/m + x1
    %      0 = mx + b => x = -b/m  = -(y1 - m x1)/m =  -y1/m + x1 OR -(y2 - m x2)/m = x2 - y2/m
    %       x = x1 - y1 (x2 - x1)/(y2 - y1);   OR x2 - y2 (x2 - x1)/(y2 - y1)
    %s(use_secant)	 = b(use_secant) - err2(use_secant) .* (b(use_secant) - a(use_secant)) ./ (err2(use_secant) - err1(use_secant));
    if any(use_secant)
        % NOTE: We only take a secant when a = c, 
        %   so it doesn't matter whether we use err2/err1 or err2/err3
        % p /q = ((a - b) * err2/err1) / ((err1 - err2)/err1)
        %      = err2 (a - b) / (err1 - err2)  -- compare to secant formula
        p1 = ab_gap .* r2;
        p(use_secant) = p1(use_secant);
        q(use_secant) = 1 - r2(use_secant);
        if verbose, fprintf('secant '), end
    end
    
    if any(try_interp)
        %pq1 = p ./ q;    % divide now (before subsetting, since it's faster); 
        %xstep(try_interp) = pq1(try_interp);  % though only for the interpolated points, so far.
        % note: this does not need subsetting because all(p(~try_interp)) = 0; 
        %  (it should be impossible for q=0 since that would imply err1 = err2, but then err isn't out of tolerance
        xstep = p ./ q;    % divide now (before subsetting, since it's faster); 
        xstep(~try_interp) = 0; % to clear nan's
    end
    
    %-----
    % 3. Override any of the above with bisection , depending on values of s
    % 3a: s is NOT between (3 * a + b)/4 and b
    %	 Brent says: new point should be up to 3/4 way between b and a
    %bi_test1a = sign(s - a) == sign(s - b);  % s is outside (a ... b); i.e. just a safety-check (and Dekker's test)
    bi_test1 = ( abs(p) >= 0.75 .* abs(ab_gap .* q) - 0.5*abs(tolx.*q) ) ; % i.e, toler

    % second test done previously
    % third test: |the current xstep| > 1/2 |xstep[n-2]|, i.e. we're moving too far away from the best. but here it's different
    %bi_test3 = abs(xstep) >= 0.5 * abs(xstep2); % or the "safer":
    bi_test3 = abs(p) >= 0.5 * abs(xstep2 .* q); % need to think about p and q here
    
    % update the use_bisection flag with the (p/q)-dependent criteria
    use_bisection = ( use_bisection | bi_test1 | bi_test3)  & err_outside_tol; % 
        
    % now accept any qualifying interpolated point:
    
    if any(use_bisection)  % it's a pretty rare event
       % s(use_bisection) = (a(use_bisection) + b(use_bisection))*0.5;
        m = -ab_gap./(2 + accel_bi); %(1+ 2.^conseqs);
        %s(use_bisection) =  b(use_bisection) - m;  % b + (a-b)/2 = (b + a)/2
        %  set xstep1 so it makes its way into xstep2 (as per Brent 1971)
        [xstep(use_bisection), xstep1(use_bisection)] =  deal(m(use_bisection)); 
        if verbose, fprintf('bisect '), end
    end

    % xstep (d in Brent) is fully updated, now compute the new guess (note xstep=0 if err is within tol)
    s = b - xstep; 
    xstep_too_small = abs(xstep) < tolx & err_outside_tol; % prevent underflow, but it converges faster if we test against tol rather than eps
    if any(xstep_too_small)
        s2 = b + sign(ab_gap) .* tolx;
        s(xstep_too_small) = s2(xstep_too_small);
    end
    %s(err_outside_tol) = s1(err_outside_tol); % preserve values that have already converged
    
    % Quick sanity check
    if ~all(use_secant | use_quad | use_bisection |  ~err_outside_tol)
        error('Somehow, we didn''t update idx: %d\n', find(~(use_secant | use_quad | use_bisection |  ~err_outside_tol)));
    end
    
    %-----
    % compute error-value for s (how far from the fixed point)
    %  this values should always be "better" than either a or b (but not necessarily better than both)
    if subsetLimit <= 0 || mean(err_outside_tol(:)) > subsetLimit
         err_s = func(s);
    else
        err_s(err_outside_tol) = func(s(err_outside_tol));
    end
    if track_fcount
        fcounter = fcounter + err_outside_tol;
    end
    
    
    %fprintf('Laggards: %d\n', sum(err_outside_tol));
    counter = counter + 1;
    if (counter > iter_limit)
        break
        error('iteration limit exceeded');
    end
    %-----------------
    %  Now reorganize a, b, c so b is in best
    %    Also: set conseqs, err_increased, ab_gap for next round
    if all(abs(err_s) < tolFn)
        % converged in y; s hold the full answer; no need to sort, etc.
        b = s;
        err2 = err_s;
        err_outside_tol=false;  % we're done!
        recompute_b = false;
    else
        % first, test that our new guess was an improvement
        %  if the new guess is larger than the old guess then most likely b is close to zero
        %  and we're just whittling away at a. (Alternatively, the function is not monotonic on [a b])
        %  Either way, we will try bisecting on the next round.
        best_is_unchanged = abs(err_s) > abs(err2); % strictly-greater performs what we really mean and in one extreme case (Heaviside) saves us from incorrect behavior
        if accelBisection
            % Furthermore, if we're just chipping away at the far side, let "bisection" divide by increasing integers.
            accel_bi = accel_bi + best_is_unchanged; 
            accel_bi(~best_is_unchanged | (accel_bi >= accelLimit)) = 0; % reset anything that behaved properly; or after 3 tries
            %      limit the number of consecutive forced-bisections (3 may be best for well-behaved fns (but may prevent convergence in singularities);
            %                                                        Inf is probably safest, though with the current method, it's never needed)
            best_is_unchanged =  accel_bi > 0; %best_is_unchanged & ; 
        end
        
        % update previous  states: (Brent's, a, fa)
        % first store prev. round's best in c: (prev round's result is in s right now)
%        d= c; err4 = err3;
        c = b;   err3 = err2; 

        %swap a,b as needed so that all b matches the sign of s
%         s_a_sign_match = (sign(err_s) == sign(err1) ) & err_outside_tol;
%         if any(s_a_sign_match)
%             % s belong in a, so move b into its place. but if we moved b in, b is now in a and c, 
%             % so lets swap 'a' into 'c' to retain three points
%             c(s_a_sign_match) = a(s_a_sign_match);      err3(s_a_sign_match)= err1(s_a_sign_match);
%             a(s_a_sign_match) = b(s_a_sign_match);      err1(s_a_sign_match)= err2(s_a_sign_match);
%             xstep1(s_a_sign_match) = xstep(s_a_sign_match); % makes a very tiny improvement in some cases
%         end
% 
%         % and now we can just copy all of s into b
%         b = s; err2 = err_s;
% 
%         % finally swap a, b if necessary, so that b is the best answer
%         err1_is_best = (abs(err2) > abs(err1)) & err_outside_tol;
%         if any(err1_is_best)
%             [a(err1_is_best), b(err1_is_best)] = deal(b(err1_is_best), a(err1_is_best));
%             [err1(err1_is_best), err2(err1_is_best)] = deal(err2(err1_is_best), err1(err1_is_best));
%         % Adding this makes it much closer to Brent's (i.e. erasing the improvement)
%         %    err_increased = ( abs(err2) >= abs(err1) );
%         end
%         
%         % very slight improvement, though makes Stage's hyperbolic case (singularity at 0) worse
%         d_closer_than_c = abs(c - b) > abs(d - b) & err_outside_tol;
%         if (d_closer_than_c)
%             c(d_closer_than_c) = d(d_closer_than_c);      err3(d_closer_than_c)= err4(d_closer_than_c);
%         end

        s_b_sign_match = (sign(err_s) == sign(err2) );
        err_s_is_best = ( abs(err_s) <= abs(err2) ) & err_outside_tol;
        a_into_b = (s_b_sign_match & ~err_s_is_best) & err_outside_tol;
        if any(a_into_b)
            % move a into b, because we're going to move s into a (note, b has already been moved into c
            b(a_into_b) = a(a_into_b);      err2(a_into_b)= err1(a_into_b);
        end

        b_into_a = (~s_b_sign_match & err_s_is_best); %  & err_outside_tol is redundant
        if any(b_into_a)
            % move b into a because we're going to move s into b; we need to save a into c (since prev-b isn't being lost, but a will)
            c(b_into_a) = a(b_into_a);      err3(b_into_a)= err1(b_into_a);
            a(b_into_a) = b(b_into_a);      err1(b_into_a)= err2(b_into_a);
        end
        % now copy s into a or b
        if any(err_s_is_best)
            b(err_s_is_best) = s(err_s_is_best);      err2(err_s_is_best)= err_s(err_s_is_best);
        end
        err_s_not_best = ~err_s_is_best & err_outside_tol;
        if any(err_s_not_best)
            a(err_s_not_best) = s(err_s_not_best);      err1(err_s_not_best)= err_s(err_s_not_best);
            xstep1(err_s_not_best) = xstep(err_s_not_best); % As per Brent; makes a very tiny improvement in some cases (& a bit worse in others?)
        end
        
        %-----
        % 0. Test if it's time to exit
        ab_gap =  (a - b); % i.e. 2 * m in Brent ( m = 0.5*(a - b) translating his paralance )
        %clear xstep so it can be partially filled w/o carryover of previous step

        %tol = eps(b) + 0.5*toler;  % BRENT: 2 * eps * abs(b)+ toler
        tolx =  2*max(1,abs(b)).*tolx_1;  % w/o max(abs(b),1.0) it fails on test f2_2
        err_outside_tol =  (0.5 .* abs(ab_gap) > tolx   &   abs(err2) > tolFn) ;
        recompute_b = true;
    end
    
    if verbose
        fprintf('%d: a: %9.5g (%9.3g), b:  %9.5g (%9.3g),c: %9.5g (%9.3g), s: %9.3g\n', fcounter, a, err1, b, err2, c, err3, err_s); % a, b, c, s);
    end
end
if recompute_b
    % should never be needed if TolX = 0 and TolFn > 0 (and at least one iter)
    err2 = func(b);
    if track_fcount
        fcounter = fcounter + 1;
    end
end
if verbose, fprintf('Brent p/q preserving C. iterations: %d; fcalls: %d; xval: %0.5g\n', counter, fcounter, b), end

// src/supporting/latin_hypercube_input.m

function latin_hypercube_input(tab, n_spectra, outdir)
    if nargin == 0
        outdir = fullfile('input', 'dataset for_verification');
        tab = readtable(fullfile(outdir, 'input_borders.csv'));
        n_spectra = 30;
    end
    out_file = fullfile(outdir, 'lh_ts.csv');
    assert(exist(out_file, 'file') == 0, '`%s` file already exists, delete it first', out_file)
    
    include = logical(tab.include);
    lb = tab.lower(include)';
    ub = tab.upper(include)';
    varnames = tab.variable(include)';

    % one row - one set of parameters
    lh = lhsdesign(n_spectra, sum(include));

    params = (ub-lb) .* lh + lb;

    if any(strcmp('LIDFa' , varnames))
        % abs(LIDFa + LIDFb) <= 1
        i_lidfa = strcmp('LIDFa', varnames);
        i_lidfb = strcmp('LIDFb', varnames);
        lidfa = params(:, i_lidfa);
        lidfb = params(:, i_lidfb);
        params(:, i_lidfa) = (lidfa + lidfb) / 2;
        params(:, i_lidfb) = (lidfa - lidfb) / 2;
    end

    t = array2table(params);
    t.Properties.VariableNames = varnames;
    t.t = (1:size(t, 1))';
    writetable(t, out_file)
    
    if verLessThan('matlab', '9.1')  % < 2016b
        varnames_in = '';
    else
        varnames_in = strjoin(varnames, ', ');
    end
    fprintf('Sampled %i parameters: %s\n', length(varnames), varnames_in)
    fprintf('Saved lut input (parameters) in `%s`\n', out_file)
end

// src/supporting/leafangles.m

function [lidf]=  leafangles(a,b)                                     
% Subroutine FluorSail_dladgen
% Version 2.3 
% For more information look to page 128 of "theory of radiative transfer models applied in optical remote sensing of
% vegetation canopies"
%
% FluorSail for Matlab
% FluorSail is created by Wout Verhoef, 
% National Aerospace Laboratory (NLR)
% Present e-mail: w.verhoef@utwente.nl
%
% This code was created by Joris Timmermans, 
% International institute for Geo-Information Science and Earth Observation. (ITC)
% Email: j.timmermans@utwente.nl
%
%% main function
F           =   zeros(1,13);
for i=1:8                                                               
    theta   =   i*10;                  %                theta_l =  10:80
    F(i)    =   dcum(a,b,theta);      %                 F(theta)
end

for i=9:12                                                              
    theta   =   80 + (i-8)*2;                         % theta_l = 82:88
    F(i)    =   dcum(a,b,theta);                     %  F(theta)
end

for i=13:13                                          %  theta_l = 90:90
    F(i) =   1;                                      %  F(theta)
end

lidf        =   zeros(13,1);
for i=13:-1:2                                                           
    lidf(i) =   F(i) -   F(i-1);                     %  lidf   =   dF/dtheta;
end
lidf(1) =   F(1);                                    %  Boundary condition

%% SubRoutines
function [F]   =  dcum(a,b,theta)
rd  =   pi/180;                                     %   Geometrical constant
if a>1 
    F    =   1 - cos(theta*rd);
else
    eps     =   1e-8;
    delx    =   1;
    
    x       =   2*rd *theta;
    theta2  =   x;
                                                                        %    
    while max(delx > eps)
        y   =   a*sin(x) + 0.5*b*sin(2*x);
        dx  =   0.5*(y - x + theta2);
        x   =   x + dx;
        delx=   abs(dx);
    end
    F    =   (2*y + theta2)/pi;                     %   Cumulative leaf inclination density function
    %pag 139 thesis says: F = 2*(y+p)/pi. 
    %Since theta2=theta*2 (in rad), this makes F=(2*y + 2*theta)/pi    
end

// src/supporting/meanleaf.m

function Fout = meanleaf(canopy,F,choice,Ps)

nl    = canopy.nlayers;
nli   = canopy.nlincl;
nlazi = canopy.nlazi;
lidf  = canopy.lidf;

% author: Dr. ir. Christiaan van der Tol (tol@itc.nl)
% date:     7   December 2007
% update:   11  February 2008 made modular (Joris Timmermans)

% update:   25 Feb 2013 Wout Verhoef : Propose name change, remove globals
%                                      and use canopy-structure for input
%
% function [F2,F3] = F1tot(F,choice,Ps)
% calculates the layer average and the canopy average of leaf properties
% per layer, per leaf angle and per leaf azimuth (36)
%
% Input:
%   F       input matrix (3D)   [nli, nlazi,nl]
%   choice  integration method  'angles'            : integration over leaf angles
%                               'angles_and_layers' : integration over leaf layers and leaf angles
%   Ps      fraction sunlit per layer [nl]
%
% Output:
%   Fout    in case of choice = 'angles': [nl]
%           in case of choice = 'angles_and_layers': [1]
Fout = zeros(nli, nlazi,nl);
switch choice
%% integration over leaf angles
    case 'angles'
        
        for j = 1:nli                                   
             Fout(j,:,:)    = F(j,:,:)*lidf(j);         % [nli, nlazi,nl]
        end
        Fout                = sum(sum(Fout))/nlazi;     % [1,1,nl]
        Fout                = permute(Fout,[3 1 2]);    % [nl]
        
%% integration over layers only         
    case 'layers'
        Fout = Ps'*F/nl;
        
%% integration over both leaf angles and layers       
    case 'angles_and_layers'
        for j = 1:nli
            Fout(j,:,:)      = F(j,:,:)*lidf(j);
        end
        
        for j = 1:nl
             Fout(:,:,j)    = Fout(:,:,j)*Ps(j);
        end
        Fout                = sum(sum(sum(Fout)))/nlazi/nl;
end

// src/supporting/plot_output_one2one.m

path_x = fullfile('output', 'SCOPE_sparse_2020-06-07-1346', 'fluxes.csv');
path_y = fullfile('output', 'sunshade_2020-06-07-2234', 'fluxes.csv');
path_y = fullfile('output', 'bigleaf_2020-06-08-0022', 'fluxes.csv');
path_y = fullfile('output', 'bigleaf_my_2020-06-08-0030', 'fluxes.csv');

opts = detectImportOptions(path_x);
df_x = readtable(path_x, opts);
df_y = readtable(path_y, opts);


flu_names = {'Rnctot','lEctot','Hctot','Actot','Tcave', ...
    'Rnstot','lEstot','Hstot','Gtot','Tsave',...
    'Rntot','lEtot','Htot','rss'};

figure
for i=1:length(flu_names)
    subplot(3, 5, i)
    flux = flu_names{i};
    x = df_x.(flux);
    y = df_y.(flux);
    plot(x, y, 'o', 'MarkerFaceColor', 'r')
%     title(flux)
    hold on
    
    % liner model
    i_nans = isnan(x);
    lm_full = fitlm(x(~i_nans), y(~i_nans));   
    lm = polyfit(x(~i_nans), y(~i_nans), 1);
    predict_x = [min(x), max(x)];
    fit = polyval(lm, predict_x);
    plot(predict_x, fit, 'r:', 'LineWidth', 2.5)  % refline(lm(1), lm(2)) % lsline()
    % metrics
    rmse = sqrt(nanmean((x - y) .^ 2));
    bias = nanmean(x - y);
    title(sprintf('%s\n%.2g (rmse), %.2g (bias), \\color{red}%.1g (r^2adj)',...
        flux, rmse, bias, lm_full.Rsquared.Adjusted))
    xlabel('SCOPE')
    ylabel('bigleaf my')
    refline(1, 0)
end

sgtitle('Bigleaf ebal no Fc vs original SCOPE')

// src/supporting/satvap.m

function es = satvap(T)

%% function [es,s]= satvap(T) 
% Author: Dr. ir. Christiaan van der Tol
% Date: 2003
%
% calculates the saturated vapour pressure at 
% temperature T (degrees C)
% and the derivative of es to temperature s (kPa/C)
% the output is in mbar or hPa. The approximation formula that is used is:
% es(T) = es(0)*10^(aT/(b+T));
% where es(0) = 6.107 mb, a = 7.5 and b = 237.3 degrees C
% and s(T) = es(T)*ln(10)*a*b/(b+T)^2

%% constants
a           = 7.5;
b           = 237.3;         %degrees C

%% calculations
es          = 6.107*10.^(7.5.*T./(b+T));
%s           = 0;%es*log(10)*a*b./(b+T).^2;

// src/supporting/slope_satvap.m

function [es,s] = slope_satvap(T)

%% function [es,s]= satvap(T) 
% Author: Dr. ir. Christiaan van der Tol
% Date: 2003
%
% calculates the saturated vapour pressure at 
% temperature T (degrees C)
% and the derivative of es to temperature s (kPa/C)
% the output is in mbar or hPa. The approximation formula that is used is:
% es(T) = es(0)*10^(aT/(b+T));
% where es(0) = 6.107 mb, a = 7.5 and b = 237.3 degrees C
% and s(T) = es(T)*ln(10)*a*b/(b+T)^2

%% constants
a           = 7.5;
b           = 237.3;         %degrees C
log10       = 2.3026;

%% calculations
es          = 6.107*10.^(7.5.*T./(b+T));
s           = es*log10*a*b./(b+T).^2;

// src/supporting/soil_respiration.m

function [R] = soil_respiration(Ts)
R = 0 + 0*Ts;    %umol m-2 s-1

// src/supporting/soltir_tp7.m

function [wn,wl,Tall] = soltir_tp7(filename)

% soltir_tp7 Reads MODTRAN tp7 file and applies a new MIT algorithm to derive
% 18 spectral functions for atmospheric correction and simulations at BOA
% and TOA

% This function reads a MODTRAN tp7 output file that contains data for 4
% runs with different albedos a:
%
% For a = 0.5 at ground altitude
% For a = 1.0 at ground altitude
% For a = 0.0 at TOA
% For a = 1.0 at TOA
%
% From the MODTRAN outputs 18 atmospheric properties are derived, 
% named t1...t18
% The radiance data are converted into units of mW m-2 sr-1 nm-1
% 
% (c) Wouter Verhoef 2011-2014

fid     = fopen(filename,'r');
modname = filename(1:size(filename,2)-4);
outname = [modname '.atm'];

for i = 1:7, fgetl(fid); end

rline = str2num(fgetl(fid));
tts   = rline(4);
cts   = cosd(tts);

s     = fgetl(fid);
rline = sscanf(s,'%10f',3);
wns   = rline(1); wne = rline(2); wstep = rline(3);
nw    = int32((wne-wns)/wstep)+1;

datarec = zeros(nw,15,4);        % MODTRAN5.2.1 tp7 15 column output format
Tall    = zeros(nw,18);          % 18 output spectra

fgetl(fid); fgetl(fid);

for ipass = 1:4
    for il=1:nw
        s=fgetl(fid);
        dline=str2num(s);
        datarec(il,:,ipass)=dline;
    end
    for j = 1:12, fgetl(fid); end
end

wn  = datarec(:,1,1);
fac = wn.*wn;
wl  = 1e7./wn;
wls = wl(nw);
wle = wl(1);

% MIT algorithm supporting fluorescence retrievals
% Wout Verhoef Sept. 2011

% Support all applications, T-18 system
% OpT in heavy absorption bands now estimnated from Planck Tb at 6500 nm
% Wout Verhoef Oct. - Nov. 2012

% Constants of Planck function

c1 = 1.191066e-22;
c2 = 14388.33;

tran_boa    = datarec(:,2,1);
tran_toa    = datarec(:,2,3);

too         = tran_toa;

toasun      = datarec(:,14,4).*fac/pi*cts;  % = Eso cos(tts) / pi

%BOA

grfl50_boa  = datarec(:,7,1).*fac;
sfem50      = datarec(:,4,1).*fac;
sfem0_boa   = 2*sfem50;
grfl100_boa = datarec(:,7,2).*fac;
delgtot_boa = grfl100_boa-sfem0_boa;
crdd        = (grfl50_boa-sfem50)./(grfl100_boa-grfl50_boa-sfem50);

rdd         = max(0,1-crdd);

OpT         = crdd.*delgtot_boa./tran_boa;

% OpT at 6500 nm is used to estimate brightness temperature of La(b), which
% is used to get a minimum radiance where O is zero

wlp         = 6500;
Lp          = interp1(wl,OpT,wlp,'nearest');
Tp          = c2/(wlp*1e-3*log(1+c1*(wlp*1e-9)^(-5)/Lp));
Lmin        = c1*(wl*1e-9).^(-5)./(exp(c2./(wl*1e-3*Tp))-1);

%TOA

grfl100_toa = datarec(:,7,4).*fac;
sfem0       = datarec(:,4,3).*fac;
delgtot_toa = grfl100_toa - sfem0;
OpTtran     = crdd.*delgtot_toa;
path100     = datarec(:,5,4).*fac;
path0       = datarec(:,5,3).*fac;

rso         = path0./toasun;

delpath     = path100 - path0;
ptem100     = datarec(:,3,4).*fac;
ptem0       = datarec(:,3,3).*fac;
delptem     = max(0,ptem100-ptem0);
delatmo     = delpath + delptem;

iT          = (wl > 4600);
ia          = (~iT & delpath == 0) | (iT & delptem == 0);

fO          = delpath./(delatmo+1e-60).*~ia + ~iT.*ia;
fT          = delptem./(delatmo+1e-60).*~ia +  iT.*ia;

O           = fO.*OpT;
T           = fT.*OpT;

% Correct cases where T = 0

i0          = (T == 0);
T(i0)       = Lmin(i0);
fT(i0)      = T(i0)./OpT(i0);
Ttran       = fT.*OpTtran;
Otran       = OpTtran-Ttran;

tdo         = delatmo./(delgtot_toa+1e-6).*tran_toa;

tdo(ia)     = 0;

gsun100_toa = datarec(:,8,4).*fac;
gsun100_boa = datarec(:,8,2).*fac;

tsstoo      = gsun100_toa./toasun;
tss         = gsun100_boa./toasun;
tsdM        = O./toasun-tss;

tsdM(ia)    = 0;

% Apply log regression to predict tsd from tdo, too, tss, and wl

Ir          = (wl>1500 & wl<1700) | (wl>2100 & wl<2300) | ...
               (wl>3950 & wl<4100);
y           = log(tsdM(Ir)./tss(Ir))-log(tdo(Ir)./too(Ir));
x           = log(wl(Ir));
n           = size(x,1);
xm          = sum(x)/n;
ym          = sum(y)/n;
a           = (x-xm)'*(y-ym)/((x-xm)'*(x-xm));
b           = ym-a*xm;
p           = a*log(wl)+b;
tsdp        = tss.*tdo./(too+1e-60).*exp(p);

% weight proportional to delatmo squared

wgtM        = delatmo.^2;

tsd         = (wgtM.*tsdM+tsdp)./(wgtM+1);

fsun        = (tss+1e-60)./(tss+1e-60+tsd);
iH          = (wl>2600);

tsdtoo      = Otran./toasun-tsstoo;
tsdtoo(iH)  = tsd(iH).*too(iH);

tssrddtoo   = tsstoo.*rdd;
toordd      = too.*rdd;
tssrdd      = tss.*rdd;
tsstdo      = fsun.*delpath.*crdd./toasun;
tsdtdo      = (1-fsun).*delpath.*crdd./toasun;
Lat         = ptem0-sfem0_boa./crdd.*tdo;
Lab         = T+sfem0_boa.*crdd;
Labtoo      = Ttran+crdd.*sfem0_boa.*too;
Labtdo      = crdd.*(delptem+sfem0_boa.*tdo);

Tall        = [toasun rso rdd tss tsd too tdo tsstoo tsdtoo tsstdo tsdtdo ...
                  tssrdd toordd tssrddtoo Lat Lab Labtoo Labtdo];

fclose(fid);

% Verification against MODTRAN

% pt0   = Lat;
% pt100 = Lat+Labtdo./(1-rdd);
% pa0   = rso.*toasun;
% pa100 = toasun.*(rso+(tsstdo+tsdtdo)./(1-rdd));
% gb50  = (toasun*.5.*(tss+(2*tsd+tssrdd)./(2-rdd))+Lab./(2-rdd)).*tran_boa;
% gb100 = (toasun.*(tss+tsd)./(1-rdd)+Lab./(1-rdd)).*tran_boa;
% gt100 = toasun.*(tsstoo+(tsdtoo+tssrddtoo)./(1-rdd))+Labtoo./(1-rdd);
% gs100 = toasun.*tsstoo;

% Write data to output file

fid=fopen(outname,'w');
str1=' WN (cm-1) WL (nm)  '; 
str2='      T1       T2       T3        T4        T5        T6        T7      ';
str3='  T8        T9        T10       T11       T12       T13       T14      ';
str4=' T15          T16          T17          T18';
str5=['                       toasun     rso      rdd       tss       tsd       too' ...
 '       tdo      tsstoo    tsdtoo    tsstdo    tsdtdo    tssrdd'  ...
 '    toordd  tssrddtoo     Lat          Lab         Labtoo       Labtdo'];
str=[str1 str2 str3 str4];
fprintf(fid,'%s\r\n',str);
fprintf(fid,'%s\r\n\r\n',str5);

for i = 1:nw
    str = sprintf('%9.2f',wn(i));
    str = [str sprintf('%10.3f',wl(i))];
    str = [str sprintf('%10.5f',Tall(i,1))];
    str = [str sprintf('%10.6f',Tall(i,2:14));];
    a   = sprintf('%14.6e',Tall(i,15:18));
    str = [str a];
    fprintf(fid,'%s\r\n',str);
end

fclose('all');

end

// src/supporting/tav.m

function Tav = tav(alfa,nr)
n2                                  =   nr.^2;
np                                  =   n2 + 1;
nm                                  =   n2 - 1;

% Stern's formula in Lekner & Dorf (1988) gives reflectance for alfa = 90 degrees

% y1 = (3*n2+2*nr+1)./(3*(nr+1).^2);
% y2 = 2*nr.^3.*(nr.^2+2*nr-1)./(np.^2.*nm);
% y3 = n2.*np.*log(nr)./nm.^2;
% y4 = n2.*nm.^2./np.^3.*log((nr.*(nr+1)./(nr-1)));

% st = y1-y2+y3-y4;

a                                   =   +((nr+1).^2)/2;
k                                   =   -((n2-1).^2)/4;
sin_a                               =   sind(alfa);
%
if alfa~=0    
    B2                              =   sin_a^2 - np/2;
    B1                              =   (alfa~=90) * sqrt( B2.^2 + k );
   
    b                               =   B1 - B2;
    b3                              =   b.^3;
    a3                              =   a.^3;
    
    ts                              =   (k.^2./(6*b3) + k./b - b./2) - ...
                                        (k.^2./(6*a3) + k./a - a./2);
                                    
    tp1                             =   -2*n2.*    (   b  -  a   ) ./ (np.^2);
    tp2                             =   -2*n2.*np.*(  log(b./a)  ) ./ (nm.^2);
    tp3                             =      n2.*    ( 1./b - 1./a ) ./ 2; 
    
%     tp4                             =   16*n2.^2.* (n2.^2+1) .* ( log(2*np.*b - nm.^2) - log(2*np.*a - nm.^2) ) ./ (np.^3.*nm.^2);    
%     tp5                             =   16*n2.^2.* (n2     ) .* ( 1./(2*np.*b - nm.^2) - 1./(2*np.*a - nm.^2)) ./ (np.^3       );

    tp4                             =	16*n2.^2.* (n2.^2+1) .* ( log((2*np.*b - nm.^2)./(2*np.*a - nm.^2))  ) ./(np.^3.*nm.^2);
    tp5                             =   16*n2.^2.* (n2     ) .* ( 1./(2*np.*b - nm.^2) - 1./(2*np.*a - nm.^2)) ./(np.^3);							 
    tp                              =   tp1 + tp2 + tp3 + tp4 + tp5;
    Tav                             =   (ts + tp) ./ (2*sin_a.^2);
else
    Tav                             =   4 *nr/((nr+1)*(nr+1));
end
return

// src/supporting/zo_and_d.m

function [zom,d] = zo_and_d(soil,canopy,constants)

% function zom_and_d calculates roughness length for momentum and zero
% plane displacement from vegetation height and LAI
%
% Date:     17 November 2008
%           17 April 2013 (structures)
%
% Author:   A. Verhoef
%           implemented into Matlab by C. van der Tol (c.vandertol@utwente.nl)
%
% Source:   Verhoef, McNaughton & Jacobs (1997), HESS 1, 81-91
% 
% usage: 
%       zo_and_d (soil,canopy)
%
% canopy fields used as inpuyt:
%   LAI         one sided leaf area index
%   hc           vegetation height (m)
%
% soil fields used:
%   Cd          Averaged drag coefficient for the vegetation              
%   CR          Drag coefficient for isolated tree
%   CSSOIL      Drag coefficient for soil
%   CD1         Fitting parameter
%   Psicor      Roughness layer correction
%
% constants used (as global)
%   kappa       Von Karman's constant
%
% output:
%   zom         roughness lenght for momentum (m)
%   d           zero plane displacement (m)
% 

%% constants
kappa   = constants.kappa;

%% parameters
CR      = canopy.CR;
CSSOIL  = soil.CSSOIL;
CD1     = canopy.CD1;
Psicor  = canopy.Psicor;
LAI     = canopy.LAI;
h       = canopy.hc;

%% calculations
sq      = sqrt(CD1*LAI/2);
G1      = max(3.3, (CSSOIL + CR*LAI/2).^(-0.5));
d       = (LAI>1E-7 & h>1E-7).*h.*(1-(1-exp(-sq))./sq);          % Eq 12 in Verhoef et al (1997)
zom     = (h-d).*exp(-kappa*G1 + Psicor);