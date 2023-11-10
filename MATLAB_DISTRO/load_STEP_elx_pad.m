function SolO = load_STEP_elx_pad(SolO,d_STEP,q_day)
    addpath(d_STEP)
    %----> Below is just to open file and dataset
    q_f   = strcat(d_STEP,q_day);
    files = dir(q_f);
    nfiles = length(files);
    disp("load_STEP_data:: I found" + num2str(nfiles) + "in this directory.")
    i =1;
    disp("load_STEP_data:: Pulling data from "+files(i).name);
    data   = spdfcdfread(files(i).name,  'ConvertEpochToDatenum',true,'Structure', true); % Read data
    %-----
    time   = data(1).Data;
    SolO.t_step   = datetime(time, 'ConvertFrom', 'datenum','InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSSSSSSSS'); % Time management (t is now in datetime)
    dss = data(51).Data; % Read average magnet channel and clean for fill values
    dss(abs(dss) >= 1e29) = NaN;
    SolO.integral_avg = dss; % Average flux from itnegral
    SolO.integral_pix = zeros([size(SolO.magnet_avg,1),size(SolO.magnet_avg,2),15]); % Initialising pixelwise fluxes froom integral
    SolO.pix_ch_txt = data(end-13).Data; % Text containing energy channel energy values
    j = 1;
    for i=6:3:50 % Spacing valid only for new data structure
        %disp(data(i).VariableName);
        dss = data(i).Data; % Read average integra; channel and clean for fill values
        dss(abs(dss) >= 1e29) = NaN; % Clean for fill however you like (recommned NaNs, fill values are \pm 1e30)
        SolO.integral_pix(:,:,j) = dss; % Load pixwelwise info (goal: Nt x Ech x pixels)
        disp("load_STEP_data:: Loading "+num2str(j) + "/ 16 pixel from ds --> " + data(i).VariableName);
        j = j+1;
    end
    % MAGNET channel
    %---> Repetition of above but for magnet channel
    dss = data(98).Data; % Read average magnet channel and clean for fill values
    dss(abs(dss) >= 1e29) = NaN;
    SolO.magnet_avg = dss;
    SolO.magnet_pix = zeros([size(SolO.magnet_avg,1),size(SolO.magnet_avg,2),15]);
    j = 1;
    for i=53:3:95
        dss = data(i).Data; % Read average magnet channel and clean for fill values
        dss(abs(dss) >= 1e29) = NaN;
        SolO.magnet_pix(:,:,j) = dss;
        disp("load_STEP_data:: Loading "+num2str(j) + "/ 16 pixel from ds --> " + data(i).VariableName);
        j = j+1;
    end
    %-----
    SolO.ele_multi = zeros(32,16); %Initialise magic numbers vector  (Ech x pixels)
    j=1;
    for i=105:120 
        SolO.ele_multi(:,j) = data(i).Data; % Read magic numbers into array
        disp(data(i).VariableName) 
        j=j+1;
    end
    SolO.elx_pix = SolO.magnet_pix.*0.0; % initialise pixelwise flux for electrons (must be same size as integral and magnet Nt x Ech x pixels )
    for it =1:size(SolO.magnet_pix,1) % Time loop
        for iE = 1:size(SolO.magnet_pix,2) % Energy loop
            for ipix = 1:size(SolO.magnet_pix,3) % Pixel loop
                SolO.elx_pix(it,iE, ipix) = (SolO.integral_pix(it,iE,ipix) - SolO.magnet_pix(it,iE,ipix))*SolO.ele_multi(iE,ipix); % Get electron flux
                %RISK: negative fluxes (waiting from news from EPD team)
            end
        end
    end
    disp("load_STEP_elx_pad:: Electron Pixelwise info loaded!")
    % /electrons - AVG (if you want) 
    % Average electron flux
    SolO.elx_flux = SolO.magnet_avg.*0.0;
    for it = 1:size(SolO.magnet_avg,1)
        for iE =1:size(SolO.magnet_avg,2)
            SolO.elx_flux(it,iE) = abs(SolO.integral_avg(it,iE) - SolO.magnet_avg(it,iE))*SolO.ele_multi(iE,16);
        end
    end
    disp("Done");
    % Last thing: build energy vector (given into the STEP data)
    elx_pix = SolO.elx_pix;
    elx_flux = SolO.elx_flux;
    E_pix = SolO.Evec_STEP; % this is the energy vector Ech from STEP (maybe add a line somewhere)
    
    SolO.elx_pix_trick = zeros(size(elx_pix));
    SolO.elx_flux_trick = zeros(size(elx_flux));
    %sSolO.EPT_trick  = zeros(size(sSolO.EPT_OMNI));
    disp(size(SolO.elx_pix_trick));
    disp(size(SolO.elx_flux_trick));
    for it =1:size(SolO.t_step,1)
        for ipix = 1:15
            SolO.elx_pix_trick(it,:,ipix) = squeeze(elx_pix(it,:,ipix)).*(E_pix.^2)'; % Mainly done for plotting purposes (Eflux \cdot E^2)
        %sSolO.EPT_trick(it,:) = squeeze(sSolO.EPT_sun(it,:)).*(sSolO.Evec_EPT.^2)';
        end
        SolO.elx_flux_trick(it,:) = squeeze(elx_flux(it,:)).*(squeeze(E_pix').^2); % Trick on the average too
    end

    % Load pixel RTN vectors (pixels x 3)
    SolO.pixel_RTN = data(end-7).Data;

    disp("Done");

end