function out = oppchanmodel(in,figs)
% Matlab implementation of the computational opponent-channel model described in:
% 
% Age-related deterioration of the representation of space in human auditory cortex
% Neurobiology of Aging
%
% Paul M. Briley (1*) & A. Quentin Summerfield (1,2)
% 1 Department of Psychology, University of York, UK
% 2 Hull York Medical School, University of York, UK
% * email: paul.briley@york.ac.uk
%
% out = oppchanmodel(in,figs);
%
% figs = 0 (no figures), 1 (channel figures, including maa predictions), 2 (ac figures), 3 (all
% figures)
%
% in = '' (default parameters), 'young' or 'youngold' or 'oldold'
% (parameters fitted to the EEG data from one of our three participant
% groups)
%
% alternatively, in can be a structure with fields chans, ac and pred
% (if a field is missing, default values are used). these fields specify
% the model parameters, see defaults for more information
%
% out is a structure containing the updated chans (channel parameters, tuning
% curves and gradients of tuning curves), ac (auditory cortex parameters and EEG response
% curves) and pred (maa predictions) fields
%

if nargin<2; figs = 3; end
if nargin<1; in = ''; end
plotacsep = 0; % if 0, collapse across cortices for plotting. otherwise, plot responses for each cortex separately

%%%% defaults (parameters fitted to EEG data from our young listener group)
chans{1}.loc = -90; % azimuth of peak response (degrees)
chans{1}.amp = 1; % peak amplitude (0-1), used for MAA predictions only. it is 1 for the young adults. for the older adults, it is the mean size of the EEG scale factor as a proportion of that for young adults
chans{1}.wid = 82; % width (degrees), corresponds to pdf standard deviation when shape is 2
chans{1}.shp = 2.6; % shape (2 = normal distribution, <2 = steeper/narrower, >2 = flatter/broader)

chans{2}.loc = 90;
chans{2}.amp = 1;
chans{2}.wid = 82;
chans{2}.shp = 2.6;

ac{1}.label = 'Left AC';
ac{1}.comp = 3; % compression factor
ac{1}.weights = [0.8907 1.1093]; % left channel, right channel, should average 1
ac{1}.scale = 88.9365; % scale factor (nAm)
ac{1}.noise = 3.2907; % noise factor (nAm)

ac{2}.label = 'Right AC';
ac{2}.comp = 3;
ac{2}.weights = [0.9934 1.0066];
ac{2}.scale = 81.4493;
ac{2}.noise = 3.3582;

pred.data = [0 5.7903 3.7489; 45 8.3496 2.7205; 60 11.3807 6.5036; 75 22.1313 12.4812]; % psychophysical data for plotting, in reference-MAA-confidenceinterval triplets
pred.k = 0.091; % conversion factor from tuning-curve gradients to MAA predictions
pred.maak = [0 5.7903]; % alternatively, give psychophysical data (reference-MAA pair) for computing conversion factor (this variable is only used if pred.k = NaN)
%%%%

%%%% plotting
chans_cols = [0 0 0; 166 166 166; 200 200 200]; chans_cols = chans_cols./255; % colours for channel tuning curves
chans_style = {'-','--','-',':'}; % line styles for channel tuning curves

sumgrad_col = [255 0 0]; sumgrad_col = sumgrad_col./255; % colour for the summed gradient across channels
sumgrad_style = '-'; % line style for summed gradient

maapred_col = [0 0 0]; maapred_col = maapred_col./255; % colour for maa predictions
maapred_style = '-'; % line style for maa predictions

maadat_col = [0 0 255]; maadat_col = maadat_col./255; % colour for maa data
maadat_style = 'o'; % line style for maa data

ac_cols = [231 120 23; 0 146 63; 0 147 221; 219 33 76; 40 22 111]; ac_cols = ac_cols./255; % colours for location-shift response curves
ac_style = {'-s','-^','-o','-+','-d'}; % line styles for response curves
ac_labs = {'-60°','-30°','0°','+30°','+60°'}; % legend labels
%%%%

ltype = 0;
if ischar(in) % simulating a listener type
    switch in
        case 'young'; ltype = 1;            
        case 'youngold'; ltype = 2;
        case 'oldold'; ltype = 3;
        case '';
        otherwise; error('unknown listener type (%s)! in should be young, youngold or oldold, or it can be a structure (see beginning of file)',in);
    end
end
if ltype
    in = loaddata(ltype);
end

if isfield(in,'chans')
    clear chans; chans = in.chans;
end
if isfield(in,'ac')
    clear ac; ac = in.ac;
end
if isfield(in,'pred')
    clear pred; pred = in.pred;
end

x = -90:1:90; % evaluate model over this range of azimuths
locs = -60:30:60; % stimuli presented from these azimuths
nlocs = length(locs);

for i = 1:length(chans) % calculate channel tuning curves and gradients
    temp = gnormpdf(x,chans{i}.loc,chans{i}.wid,chans{i}.shp); % get tuning curve    
    temp = temp - min(temp); % these lines make the curve extend from 0 to chans{i}.amp
    temp = temp./max(temp);
    chans{i}.tun = temp.*chans{i}.amp;
    clear temp;
    
    temp1 = [diff(chans{i}.tun) 0]; % calculate gradient of channel tuning curve
    temp2 = -fliplr([diff(fliplr(chans{i}.tun)) 0]);
    temp = abs(temp1) + abs(temp2);
    temp(2:end-1) = temp(2:end-1)./2;
    chans{i}.grd = temp;
    clear temp; clear temp1; clear temp2;
    
    if i==1; grads = chans{i}.grd; else grads = grads + chans{i}.grd; end % summing gradients across channels
end

for i = 1:length(ac); ac{i}.resp = zeros(nlocs); end % pre-shift location x post-shift location
for pre_loc = 1:nlocs % calculate response sizes for each combination of pre- and post-shift location
    for post_loc = 1:nlocs
        for i = 1:length(chans)
            temp = chans{i}.tun;
            temp = temp./max(temp); % normalise channel to peak at 1 (the channel amplitude applies only to the psychophysical predictions; for the EEG data, the EEG scale factors control overall response size)
            pre_resp = temp(x==locs(pre_loc));
            post_resp = temp(x==locs(post_loc));
            d = post_resp - pre_resp;
            d(d<0) = 0;
            for ii = 1:length(ac)
                ac{ii}.resp(pre_loc,post_loc) = ac{ii}.resp(pre_loc,post_loc) + d.*ac{ii}.weights(i);
            end
        end
    end
end

for i = 1:length(ac) % apply compression, scale response sizes and add the noise factor
    if ac{i}.comp; ac{i}.resp = 1 - exp(-ac{i}.resp.*ac{i}.comp); end % an alternative form of compression is ac{i}.resp = ac{i}.resp.^ac{i}.comp (where ac{i}.comp is between 0 and 1)
    ac{i}.resp = ac{i}.resp .* ac{i}.scale + ac{i}.noise;
end

pred.refs = x; % compute maa predictions based on summed tuning-curve gradients
pred.grads = grads;
pred.igrads = 1./grads;
if isnan(pred.k)
    pred.k = pred.maak(2)./pred.igrads(x==pred.maak(1));
end
pred.maas = pred.igrads.*pred.k;

if sum(figs==[1 3]) % channel figures
    figure; % tuning curves
    for i = 1:length(chans)
        plot(x,chans{i}.tun,chans_style{i},'Color',chans_cols(i,:),'linewidth',2); hold on;
    end
    title('Channel tuning curves');
    xlabel('Azimuth in degrees');
    ylabel('Response magnitude')
    axis([x(1) x(end) 0 1]);
    plot([0 0],[0 1],'k--');
    
    figure; % tuning curve gradients
    for i = 1:length(chans)
        plot(x,chans{i}.grd,chans_style{i},'Color',chans_cols(i,:),'linewidth',2); hold on;
    end
    plot(x,grads,sumgrad_style,'Color',sumgrad_col,'linewidth',2); % summed gradient across channels
    title('Tuning curve gradients');
    xlabel('Azimuth in degrees');
    ylabel('Response magnitude')
    mx = get(gca,'YLim');
    axis([x(1) x(end) 0 mx(2)]);
    plot([0 0],[0 mx(2)],'k--');
    
    figure; % maa predictions
    plot(x,pred.maas,maapred_style,'Color',maapred_col,'linewidth',2);
    title('MAA predictions');
    xlabel('Reference azimuth in degrees');
    ylabel('MAA in degrees');
    hold on;
    errorbar(pred.data(:,1),pred.data(:,2),pred.data(:,3),maadat_style,'Color',maadat_col);
    set(gca,'XLim',[pred.data(1,1)-5 pred.data(end,1)+5]);
    %mx = get(gca,'YLim'); if mx(2)>60; mx(2) = 60; end
    mx = [0 60];
    set(gca,'YLim',[0 mx(2)]);
    plot(get(gca,'XLim'),[30 30],'k--');
end

if sum(figs==[2 3]) % ac figures
    for i = 1:length(ac) % average across cortices (for if plotacsep = 0)
        if i==1; temp = ac{i}.resp; else temp = temp + ac{i}.resp; end
    end
    temp = temp./length(ac);
    if plotacsep; nac = length(ac); else nac = 1; end
    
    for i = 1:nac % location-shift responses, using absolute pre-shift location
        figure;
        h = zeros(1,nlocs); 
        for ii = 1:nlocs % for each post-shift location
            if plotacsep; temp = ac{i}.resp; end              
            h(ii) = plot(locs,temp(:,ii),ac_style{ii},'Color',ac_cols(ii,:));
            hold on;
        end
        if plotacsep; title(ac{i}.label); else title('Averaged across cortices'); end
        leg = legend(h,ac_labs,'location','SouthEast');
        set(get(leg,'title'),'String','Post-shift location:');
        mx = get(gca,'YLim');
        axis([locs(1) locs(end) mx(1) mx(2)]);
        xlabel('Absolute pre-shift location in degrees');
        ylabel('Response magnitude in nAm');
    end
    
    for i = 1:nac % location-shift responses, using relative pre-shift location
        figure;
        h = zeros(1,nlocs);
        for ii = 1:nlocs % for each post-shift location
            if plotacsep; temp = ac{i}.resp; end
            mylocs = locs - locs(ii);
            h(ii) = plot(mylocs,temp(:,ii),ac_style{ii},'Color',ac_cols(ii,:));
            hold on;
        end
        if plotacsep; title(ac{i}.label); else title('Averaged across cortices'); end
        leg = legend(h,ac_labs,'location','SouthEast');
        set(get(leg,'title'),'String','Post-shift location:');
        mx = get(gca,'YLim');
        axis([locs(1)-locs(end) locs(end)-locs(1) mx(1) mx(2)]);
        xlabel('Relative pre-shift location in degrees');
        ylabel('Response magnitude in nAm');
    end
end
out.chans = chans;
out.ac = ac;
out.pred = pred;

function out = gnormpdf(x,m,s,b)
% m is azimuth of peak response
% s is width
% b is shape parameter
s = s.*sqrt(2); % included to make s correspond to standard deviation when b = 2 (normal distribution)
part1 = b./(2.*s.*gamma(1./b));
part2 = exp(-(abs(x-m)./s).^b);
out = part1.*part2;

function in = loaddata(ltype)
switch ltype
    case 1 % young
        chans{1}.loc = -90; chans{1}.amp = 1; chans{1}.wid = 82; chans{1}.shp = 2.6;
        chans{2}.loc = 90; chans{2}.amp = 1; chans{2}.wid = 82; chans{2}.shp = 2.6;
        ac{1}.label = 'Left AC'; ac{1}.comp = 3; ac{1}.weights = [0.8907 1.1093]; ac{1}.scale = 88.9365; ac{1}.noise = 3.2907;
        ac{2}.label = 'Right AC'; ac{2}.comp = 3; ac{2}.weights = [0.9934 1.0066]; ac{2}.scale = 81.4493; ac{2}.noise = 3.3582;
        pred.data = [0 5.7903 3.7489; 45 8.3496 2.7205; 60 11.3807 6.5036; 75 22.1313 12.4812]; pred.k = 0.091;
    case 2 % youngold
        chans{1}.loc = -90; chans{1}.amp = 0.7299; chans{1}.wid = 87; chans{1}.shp = 2.2;
        chans{2}.loc = 90; chans{2}.amp = 0.7299; chans{2}.wid = 87; chans{2}.shp = 2.2;
        ac{1}.label = 'Left AC'; ac{1}.comp = 3; ac{1}.weights = [0.8840 1.1160]; ac{1}.scale = 65.1448; ac{1}.noise = 6.3310;
        ac{2}.label = 'Right AC'; ac{2}.comp = 3; ac{2}.weights = [1.0442 0.9558]; ac{2}.scale = 59.2149; ac{2}.noise = 5.4591;
        pred.data = [0 6.1331 2.2366; 45 10.2045 3.7220; 60 9.1794 2.9366; 75 19.1181 8.0054]; pred.k = 0.091;
    case 3 % oldold
        chans{1}.loc = -90; chans{1}.amp = 0.4441; chans{1}.wid = 73; chans{1}.shp = 3.6;
        chans{2}.loc = 90; chans{2}.amp = 0.4441; chans{2}.wid = 73; chans{2}.shp = 3.6;
        ac{1}.label = 'Left AC'; ac{1}.comp = 2; ac{1}.weights = [0.9619 1.0381]; ac{1}.scale = 48.6382; ac{1}.noise = 6.7692;
        ac{2}.label = 'Right AC'; ac{2}.comp = 2; ac{2}.weights = [0.8893 1.1107]; ac{2}.scale = 27.0258; ac{2}.noise = 8.0610;
        pred.data = [0 8.3319 2.6899; 45 22.1833 6.5283; 60 42.3492 13.1188; 75 39.9189 16.1739]; pred.k = 0.091;
end
in.chans = chans;
in.ac = ac;
in.pred = pred;
