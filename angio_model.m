% Angiogenesis Model
% UCI COSMOS 2023:  Tissue and Tumor Modeling (Cluster 3)
% Provided by: Huaming Yan

function angio_model
global use_ginput tumor_r domain tips traces vgp frames N_grid_pts ... % tumor_ctr
    tip_min_dist N_tip N_tumor tumor_min_dist atm_min_dist is_in_device ...
    branch_tip_min branch_tip_max branch_tip_vegf v_rand tip_is_branch ...
    rotate_by tip_delta tips_init tip_anas_after sp device_k dx nu tm ve ...
    tumor_ctr_level tumor_D tumor_lam tumor_d tumor_dt ...
    nu_D nu_d nu_p nu_u nu_pd nc nc_v nu_dt ...
    vegf_D vegf_u vegf_d vegf_p0 vegf_p vegf_v0 vegf_dt ...
    branch_ves_min branch_ves_max branch_ves_vegf 

%% Configuration
% General parameters
use_ginput = 1;         % use mouse cursor to pick tumor and vessel seeds
rng(0);                 % random seed
frames = 1000;          % number of time steps
domain = [-1 1];        % x-axis of the device
device_k = 3/4;         % slope of device boundary
new_rand_every = 25;    % change randomness every few time steps
N_grid_pts = 1e2;       % number of grid points for continuous equations
 
% Initial condition, if not picked manually
N_tip = 10;             % number of initial tip cells
N_tumor = 5;            % number of initial tumors
tumor_r = 0.03;         % initial tumor size
tip_min_dist = 0.3;     % minimum distance between two initial tip cells
tumor_min_dist = 0.3;   % minimum distance between two initial tumors
 
% Tip cell evolution
dx = 0.005;             % distance that a tip cell travels per time step
atm_min_dist = 0.001;   % min. dist between tip cells and vessels before anastomosis
v_rand = 0.5;           % magnitude of random velocity relative to the VEGF gradient
tip_anas_after = 0.02;  % no anastomosis for a new branch before it travels this far

% Vessel branching
new_tip_freq = 0;       % a new seed every few iterations. 0 to disable
branch_tip_min = 0.005; % min probability of tip cell branching
branch_tip_max = 0.01;  % max probability
branch_tip_vegf = 0.05; % positive gain by VEGF
branch_ves_min = 0.005; % min probability of branching along vessels
branch_ves_max = 0.08;  % max probability
branch_ves_vegf = 0.1;  % positive gain by VEGF
tip_delta = 60;         % initial angle (in degrees) between two branches
 
% Tumor equation
tumor_D = 0.001;        % diffusivity
tumor_lam = 1.0;        % mitosis rate
tumor_d = 0.3;          % death rate
tumor_dt = 0.02;        % time step

% Nutrient equation
nu_D = 0.01;            % diffusivity 
nu_d = 0.0;             % decay
nu_u = 0.2;             % uptake by tumors
nu_p = 0.1;             % production by the host
nu_pd = 3.0;            % production by the vessels
nc = 0.5;               % background nutrient level
nc_v = 2.0;             % nutrient level in the vessels
nu_dt = 0.01;           % nutrient time step
 
% VEGF equation
vegf_D  = 0.01;         % diffusivity
vegf_u  = 6.0;          % uptake by vessels
vegf_d  = 0.005;        % natural decay
vegf_p0 = 0.001;        % background production
vegf_p  = 0.5;          % production by tumors
vegf_v0 = 1.0;          % maximum concentration
vegf_dt = 0.01;         % time step

% Utilities
tumor_ctr_level = 0.3;  % tumor contour level
ylims = [20 2 100];     % y-axis limit of [vessel_len, tumor_vol, num_branches]
make_movie = 0;         % make a movie. Memory expensive and slows down the evolution!

is_in_device = @(x) abs(x(1))+abs(x(2)/device_k) <= 1;
rotate_by = @(t) [cos(t) -sin(t); sin(t) cos(t)];

%% Initialization
% continuous computational domain
sp = linspace(domain(1),domain(2),N_grid_pts);
[X,Y] = meshgrid(sp,sp);

% seed tumors and tip cells
if exist('use_ginput','var') && use_ginput
    ginput_seeds;
else
    seed_tip_cells;
    seed_tumor;
end
init_cont_field;

% statistics
vessel_len = zeros(frames,1); 
tumor_vol = zeros(frames,1); tumor_vol(1) = sum(sum(tm))*(sp(2)-sp(1))^2;
num_branch = zeros(frames,1); num_branch(1) = N_tip;

%% Set up the plot
figure('position',[500 400 1300 500]);
set(gcf,'color',[1 1 1]*0.8)
subplot(2,3,1); hold on
% set(gca,'Position',[0.1 0.45 0.2 0.5])
axis off; axis square
xlim(domain); ylim(domain)
plot_device_boundary;

% plot the VEGF field
vegf = field_only_in_device(ve);
vegf = vegf/max(max(vegf));
hv = pcolor(X,Y,flipud(rot90(vegf,1)));
% caxis([0 1]);
set(hv,'FaceColor','interp'); set(hv,'EdgeColor','none');
hc = colorbar; set(get(hc,'title'),'string','VEGF','Fontsize',14); set(hc,'YTick',[])
load('fjhot5','fjhot'); colormap(gca,fjhot)

% % plot the initial tumors
% htm = cell(1,N_tumor);
% for i=1:N_tumor
%     htm{i} = rectangle('Position',[tumor_ctr(i,:)-tumor_r tumor_r*2 tumor_r*2], ...
%         'Curvature',[1 1],'EdgeColor','g','LineWidth',2);
% end

% show legends
% rectangle('Position',[domain(1)+0.02 domain(1)*0.65 tumor_r*2 tumor_r*2], ...
%     'Curvature',[1 1],'EdgeColor','g','LineWidth',2)
% text(domain(1)+0.15,domain(1)*0.62, 'Tumors', 'Fontsize',12)
plot([domain(1) domain(1)+0.1],[domain(1)*0.63 domain(1)*0.63], ...
    'r','linewidth',2)
text(domain(1)+0.15,domain(1)*0.62, 'Vessels', 'Fontsize',12)
rectangle('Position',[domain(1)+0.035 domain(1)*0.74 0.02 0.02], ...
    'Curvature',[1 1],'FaceColor',[0.95 0.95 0],'EdgeColor',[0.95 0.95 0])
text(domain(1)+0.15,domain(1)*0.72, 'Tip cell seeds', 'Fontsize',12)
if branch_tip_max>0
    rectangle('Position',[domain(1)+0.035 domain(1)*0.84 0.02 0.02], ...
        'Curvature',[1 1],'FaceColor','k','EdgeColor','k')
    text(domain(1)+0.15,domain(1)*0.82, 'Branching points', 'Fontsize',12)
end
ht = text(0.5*domain(2),0.8*domain(1),sprintf('0/%d',frames),'Fontsize',12);

% show the nutrient field
subplot(2,3,2); hold on
axis off; axis square
nuf = field_only_in_device(nu);
nuf = nuf/max(max(nuf));
hn = pcolor(X,Y,flipud(rot90(nuf,1))); 
% caxis([0 1]);
hc = colorbar; set(get(hc,'title'),'string','Nutrients','Fontsize',14); set(hc,'YTick',[]) 
plot_device_boundary;
load('fjhot','fjhot'); colormap(gca,fjhot)
set(hn,'FaceColor','interp')
set(hn,'EdgeAlpha',0)

% show the tumor field
subplot(2,3,3); hold on
% set(gca,'Position',[0.4 0.45 0.2 0.5])
axis off; axis square
tmf = field_only_in_device(tm);
if ~N_tumor
    htm = pcolor(X,Y,tmf*NaN); 
else
    htm = pcolor(X,Y,flipud(rot90(tmf,1))); 
    % caxis([0 1]);
    [~, htmf] = contour(X,Y,flipud(rot90(tmf,1)),[1 1]*tumor_ctr_level,'g','LineWidth',2);
    hc = colorbar; set(get(hc,'title'),'string','Tumor','Fontsize',14); set(hc,'YTick',[])
end
plot_device_boundary;
load('fjhot','fjhot'); colormap(gca,fjhot)
set(htm,'FaceColor','interp')
set(htm,'EdgeAlpha',0)
if N_tumor
    plot([domain(1) domain(1)+0.1],[domain(1)*0.73 domain(1)*0.73], 'g','linewidth',2)
    text(domain(1)+0.15,domain(1)*0.72, [num2str(tumor_ctr_level) ' contour'], 'Fontsize',12)
end

% show the statistics
subplot(2,3,4); hold on
set(gca,'Position',[0.1 0.1 0.15 0.2])
gcv = gca; set(gcv,'FontSize',12);
xlabel('Time'); xlim([0 frames]);
hvl = plot(1,vessel_len(1),'linewidth',2);
title('Vessel length'); ylim([0 ylims(1)])

subplot(2,3,5); hold on
set(gca,'Position',[0.4 0.1 0.15 0.2])
gcb = gca; set(gcb,'FontSize',12);
xlabel('Time'); xlim([0 frames]);
hbn = plot(1,num_branch(1),'linewidth',2);
title('Number of branches'); ylim([0 ylims(3)])

subplot(2,3,6); hold on
set(gca,'Position',[0.7 0.1 0.15 0.2])
gct = gca; set(gct,'FontSize',12);
xlabel('Time'); xlim([0 frames]);
htv = plot(1,tumor_vol(1),'linewidth',2);
title('Tumor volume'); ylim([0 ylims(2)])

% create a line object for each vessel
subplot(2,3,1); hold on
hp = cell(N_tip,1);
for i=1:N_tip
    hp{i} = plot(traces(i,1,1),traces(i,1,2),'r','linewidth',2);
end
tip_init_shown = zeros(N_tip,1);

subplot(2,3,2); set(gca,'Position',[0.35 0.35 0.25 0.55])
subplot(2,3,3); set(gca,'Position',[0.65 0.35 0.25 0.55])
subplot(2,3,1); set(gca,'Position',[0.05 0.35 0.25 0.55])

if exist('make_movie','var') && make_movie
    images = cell(frames,1); % store each frame to make a movie
end

%% Evolve the tip cells and show the traces
tc_v = zeros(N_tip,2); % tip cell velocity
tc_stopped = zeros(N_tip,1); % records if a tip cell is stopped
for time=2:frames
    % create new tip cells
    if new_tip_freq>0 && ~mod(time,new_tip_freq)
        create_new_tip_cell(time);
        
        % update the velocity and records
        tc_v(N_tip,:) = get_tip_velocity(N_tip);
        tc_stopped(N_tip) = false;
        tip_init_shown(N_tip) = false;
        vgp(N_tip,:,:) = zeros(N_grid_pts);
        hp{N_tip} = plot(traces(N_tip,time,1),traces(N_tip,time,2),'r','linewidth',2);
    end
    
    % create new branches at tip cells
    for i=1:N_tip
        if check_new_branch_at_tip(i)
            create_new_branch_at_tip(i,time,tc_v(i,:));
            
            % update the velocity of the old and new tip cell
            tc_v(N_tip,:) = rotate_by(-tip_delta/180*pi) * tc_v(i,:)';
            tc_v(i,:) = rotate_by(tip_delta/180*pi) * tc_v(i,:)';
            
            % update the records
            tc_stopped(N_tip) = false;
            tip_init_shown(N_tip) = false;
            vgp(N_tip,:,:) = zeros(N_grid_pts);
            hp{N_tip} = plot(traces(N_tip,time,1),traces(N_tip,time,2),'r','linewidth',2);
        end
    end
    
    % create new branches along vessels
    for i=1:N_tip
        [created,vd] = check_create_new_branch_along_vessel(i,time);
        if created
            tc_v(N_tip,:) = vd;
            
            % update the records
            tc_stopped(N_tip) = false;
            tip_init_shown(N_tip) = false;
            vgp(N_tip,:,:) = zeros(N_grid_pts);
            hp{N_tip} = plot(traces(N_tip,time,1),traces(N_tip,time,2),'r','linewidth',2);
        end
    end
    num_branch(time) = N_tip;
    
    % evolve tip cells
    vessel_len_incre = 0;
    for i=1:N_tip
        % if a tip cell is in a tumor, or reaches another vessel, 
        % then stop its movement
        if tc_stopped(i), continue; end
        if check_tc_stop(i), tc_stopped(i) = true; continue; end
        
        % generate new random velocity
        if ~mod(time,new_rand_every)
            tc_v(i,:) = get_tip_velocity(i);
        end

        % if the tip cell will self-intersect, stop
        if check_self_intersection(i,tc_v(i,:))
            tc_stopped(i) = true;
            continue
        end
        
        tips(i,:) = tips(i,:) + tc_v(i,:) * dx;
        vessel_len_incre = vessel_len_incre + dx; % tc_v(i,:) is a unit vector
    end
    traces(:,time,:) = tips;
    
    % update the grid points 'occupied' by vessels
    maintain_vgp;
    
    % solve the continuous equations
    update_cont_field;
    vegf = field_only_in_device(ve);
    nuf = field_only_in_device(nu);
    if N_tumor, tmf = field_only_in_device(tm); end
    set(hv,'CData',flipud(rot90(vegf,1)));
    set(hn,'CData',flipud(rot90(nuf,1)));
    if N_tumor
        set(htm,'CData',flipud(rot90(tmf,1)));
        set(htmf,'ZData',flipud(rot90(tmf,1)));
    end
    
    % plot the trace of tip cells. Change XData and YData of the handles
    for i=1:N_tip
        set(hp{i},'XData',traces(i,1:time,1),'YData',traces(i,1:time,2));
    end
    
    % show the initial location of tip cells
    nbt = 0;
    for i=1:N_tip
        if ~tip_init_shown(i) && ~isnan(traces(i,time,1))
            % Show the initial position of tip cells, if the tip cell has
            % been generated.
            if tip_is_branch(i), tc = 'k'; else, tc = 'y'; end
            rectangle('Position',[tips_init(i,:) 0.02 0.02], ...
                'Curvature',[1 1],'FaceColor',tc,'EdgeColor',tc)
            tip_init_shown(i) = true;
            nbt = nbt + 1;
        end
    end
    num_branch(time) = num_branch(time-1) + nbt*(time>2);
    
    % update the statistics
    vessel_len(time) = vessel_len(time-1) + vessel_len_incre;
    tumor_vol(time) = sum(sum(tm))*(sp(2)-sp(1))^2;
    set(hvl,'XData',1:time,'YData',vessel_len(1:time));
    set(htv,'XData',1:time,'YData',tumor_vol(1:time));
    set(hbn,'XData',1:time,'YData',num_branch(1:time));
    yl = get(gcv,'YLim'); if vessel_len(time)>yl(2), set(gcv,'YLim', yl+[0 5]); end
    yl = get(gct,'YLim'); if tumor_vol(time)>yl(2), set(gct,'YLim', yl+[0 1]); end
    yl = get(gcb,'YLim'); if num_branch(time)>yl(2), set(gcb,'YLim', yl+[0 100]); end
    
    % show the progress
    set(ht,'string',sprintf('%d/%d',time,frames));
    drawnow nocallbacks
    
    % store each frame to make a movie
    if exist('make_movie','var') && make_movie
        images{time} = getframe(gcf); 
    end
end

% save the statistics and a plot of vessels at final time
save('angio_stat','traces','vessel_len','tumor_vol','num_branch')
h = figure('visible','off'); hold on
axis off; axis square
xlim(domain); ylim(domain)
plot_device_boundary;
for i=1:N_tip
    plot(traces(i,1:frames,1),traces(i,1:frames,2),'r','linewidth',2);
    if isnan(tips_init(i,1)), continue; end
%     rectangle('Position',[tips_init(i,:) 0.02 0.02], ...
%         'Curvature',[1 1],'FaceColor','r','EdgeColor','r')
end
print(h,'-dpng','vessels','-r300');
close(h)

% make a movie of the evolution
if exist('make_movie','var') && make_movie
    disp('Making movie')
    writerObj = VideoWriter('vessels.avi');
    writerObj.FrameRate = 30;
    open(writerObj);
    for t=1:frames
        if isempty(images{t}), continue; end
        thisimage = images{t}.cdata;
        writeVideo(writerObj, thisimage);
    end
    close(writerObj);
    if exist('ffmpeg.m','file'), ffmpeg('vessels.avi'); end
end
end

function plot_device_boundary
global device_k

plot([-1 0],[0 device_k],'k')
plot([-1 0],[0 -device_k],'k')
plot([0 1],[device_k 0],'k')
plot([0 1],[-device_k 0],'k')
end

function ginput_seeds
global N_tip tips tumor_ctr N_tumor domain is_in_device

h = figure; hold on
set(gcf,'color',[1 1 1]*0.8)
axis square
axis off
plot_device_boundary;
rectangle('Position',[domain(1)+0.035 domain(1)*0.84 0.04 0.04], ...
    'Curvature',[1 1],'FaceColor',[0.95 0.95 0],'EdgeColor',[0.95 0.95 0])
text(domain(1)+0.15,domain(1)*0.82, 'Tip cell seeds', 'Fontsize',12)
rectangle('Position',[domain(1)+0.035 domain(1)*0.94 0.04 0.04], ...
    'Curvature',[1 1],'FaceColor','g','EdgeColor','g')
text(domain(1)+0.15,domain(1)*0.92, 'Tumors', 'Fontsize',12)

waitfor(msgbox('Use mouse cursor to pick tip cell seeds. Press enter to finish.'))
tips = [];
while true
    [x,y] = ginput(1);
    if isempty(x), break; end
    if ~is_in_device([x y]), disp('Out-of-device location ignored'); continue; end
    tips = [tips; x y];
    scatter(x,y,'yo','filled');
end
N_tip = size(tips,1);
initialize_tip_cells;

waitfor(msgbox('Use mouse cursor to pick tumor locations. Press enter to finish.'))
tumor_ctr = [];
while true
    [x,y] = ginput(1);
    if isempty(x), break; end
    if ~is_in_device([x y]), disp('Out-of-device location ignored'); continue; end
    tumor_ctr = [tumor_ctr; x y];
    scatter(x,y,'go','filled');
end
N_tumor = size(tumor_ctr,1);

close(h);
end

%% VEGF related
function c = vegf_conc_init(x,y)
global tumor_r tumor_ctr %is_in_device

% Note: when using this as the initial condition, do not include NaN!
% % 0.98 allows the sawtooth shape to lie on the boundary line
% if ~is_in_device([x y]*0.98), c = NaN; return; end

diff_vegf = 1;

% % Gaussian source
% c = 0;
% for i=1:size(tumor_ctr,1)
%     expo = (x-tumor_ctr(i,1))^2 + (y-tumor_ctr(i,2))^2;
%     c = c + exp(-expo/diff_vegf);
% end

% from the nevi sims
ra = tumor_r*0.5; % 50% of the tumor producing the factor
c = 0;
for iv=1:size(tumor_ctr,1)
    r = norm([x y]-tumor_ctr(iv,:));
    if r<tumor_r
        % inside a tumor
        c = c + besseli(0,tumor_r/diff_vegf) ...
            * (-tumor_r*besselk(1,tumor_r/diff_vegf) ...
            + ra*besselk(1,ra/diff_vegf));
    else
        % outside a tumor
        c = c + (tumor_r*besseli(1,tumor_r/diff_vegf) ...
            - ra*besseli(1,ra/diff_vegf))*besselk(0,r/diff_vegf);
    end
end
end

function c = vegf_conc(x,y)
global ve

ix = get_vgp([x y]);
c = ve(ix(1),ix(2));
end

function dv = vegf_gradient(x,y)
h = 1e-2;

dv = zeros(1,2);
dv(1) = (vegf_conc(x+h,y) - vegf_conc(x-h,y)) / (2*h);
dv(2) = (vegf_conc(x,y+h) - vegf_conc(x,y-h)) / (2*h);
end

% get the grid point 'occupied' by (i.e. closest to) x
function ix = get_vgp(x)
global sp

ix = floor((x-sp(1))/(sp(2)-sp(1)));
ix = min(max(ix,1),length(sp));
end

%% Tip cell related
% random walk corrected to VEGF gradient, every few time steps
function v = get_tip_velocity(i)
global tips v_rand

v = vegf_gradient(tips(i,1),tips(i,2));
rnd = randn(1,2);
v = v/norm(v) + rnd/norm(rnd)*v_rand;
v = v/norm(v);
end
            
% generate a random seeding position within the device
function tp = gen_seed_pos
global domain is_in_device

while true
    tp = domain(1) + (domain(2)-domain(1)) * rand(1,2); % proposed position
    % 0.95 ensures some distance from the boundary
    if is_in_device(tp/0.95), return; end
end
end

function seed_tip_cells
global N_tip domain tips tip_min_dist 

tips = repmat(domain*2,[N_tip 1]); % initial out-of-domain positions
for i=1:N_tip
    k = 1;
    while k<1e4
        tp_rejected = false;
        tp = gen_seed_pos;
        for j=1:N_tip
            if i==j, continue; end
            if norm(tips(j,:)-tp)<tip_min_dist
                tp_rejected = true;
                break
            end
        end
        if ~tp_rejected, break; end
    end
    if tp_rejected, disp('seeding tip cells failed'); return; end
    tips(i,:) = tp;
end

initialize_tip_cells;
end

function initialize_tip_cells
global N_tip tips tips_init traces frames tip_is_branch vgp N_grid_pts

% Initial locations of tip cells
tips_init = tips;

% Traces of tip cells. Each row is the trace of a tip cell
traces = ones(N_tip, frames, 2)*NaN;
traces(:,1,:) = tips;

% Whether a tip cell is a branching point
tip_is_branch = zeros(N_tip,1);

% Grid points 'occupied' by (i.e. closest to) vessels
vgp = zeros(N_tip,N_grid_pts,N_grid_pts);
maintain_vgp;
end

function create_new_tip_cell(time)
global N_tip tips tips_init traces frames tip_min_dist tip_is_branch

k = 1;
while k<1e4
    tp_rejected = false;
    tp = gen_seed_pos;
    for i=1:N_tip
        if norm(tips(i,:)-tp)<tip_min_dist
            tp_rejected = true;
            break
        end
    end
    if ~tp_rejected, break; end
end
if tp_rejected, disp('seeding tip cells failed'); return; end

N_tip = N_tip + 1;
tips(N_tip,:) = tp;
tips_init(N_tip,:) = tp;
traces(N_tip,:,:) = ones(frames, 2)*NaN;
traces(N_tip,time,:) = tp;
tip_is_branch(N_tip) = false;
end

function create_new_branch_at_tip(id,time,v)
global N_tip tips tips_init traces frames tip_is_branch atm_min_dist ...
    tip_delta rotate_by

N_tip = N_tip + 1;

% create new tip cell at a distance away from the old cell
tips(N_tip,:) = tips(id,:) + atm_min_dist*1.5 ...
    * (rotate_by(-tip_delta/180*pi) * v')';

tips_init(N_tip,:) = tips(N_tip,:);
traces(N_tip,:,:) = ones(frames, 2)*NaN;
traces(N_tip,time,:) = tips(N_tip,:);
tip_is_branch(N_tip) = true;
end

function [created,vd] = check_create_new_branch_along_vessel(id,time)
global N_tip tips tips_init traces frames tip_is_branch atm_min_dist rotate_by ...
    branch_ves_min branch_ves_max branch_ves_vegf 

if branch_ves_max<=0, created = false; vd = []; return; end

% % if the tip cell has been stopped, do not branch
% DISABLED. the vessel can still branch
% if check_tc_stop(id), created = false; vd = []; return; end

% pick a random point along the vessel
tm = randi(time);
loc = [traces(id,tm,1) traces(id,tm,2)];
if isnan(loc(1)), created = false; vd = []; return; end

% check branching probability
vegf = vegf_conc(loc(1),loc(2));
prob = branch_ves_min + (branch_ves_max-branch_ves_min) ...
    * branch_ves_vegf*vegf / (1+branch_ves_vegf*vegf);
if rand>prob, created = false; vd = []; return; end

N_tip = N_tip + 1;

% create new tip cell at a distance away from the vessel
if tm>1 && tm<frames && ~isnan(traces(id,tm+1,1))
    % branch is perpendicular to the existing vessel
    vd = [traces(id,tm+1,1) traces(id,tm+1,2)] - loc;
    vd = (rotate_by(-pi/2) * vd')'/norm(vd);
else
    vd = rand(1,2);
    vd = vd/norm(vd);
end
tips(N_tip,:) = loc + atm_min_dist*1.5 * vd;
    
tips_init(N_tip,:) = tips(N_tip,:);
traces(N_tip,:,:) = ones(frames, 2)*NaN;
traces(N_tip,time,:) = tips(N_tip,:);
tip_is_branch(N_tip) = true;

created = true;
end

function stop = check_tc_stop(id)
global N_tip tips vgp ... %traces atm_min_dist frames ... % tumor_ctr tumor_r
    is_in_device tip_is_branch tips_init tip_anas_after %tm tumor_ctr_level
stop = false;

% if tip cell is out of device, stop
if ~is_in_device(tips(id,:))
    stop = true;
    return
end

% % if tip cell is in a tumor, stop
% dist_to_tumor = vecnorm((tumor_ctr-tips(id,:))');
% if sum(dist_to_tumor<tumor_r)
%     stop = true;
%     return
% end

% % if tip cell is in a tumor, stop (continuous version)
% ix = get_vgp(tips(id,:));
% if tm(ix(1),ix(2))>tumor_ctr_level
%     stop = true;
%     return
% end

% if tip cell crosses the path of another vessel, stop
% unless it's a branch that hasn't traveled far
if tip_is_branch(id) && norm(tips(id,:)-tips_init(id,:))<tip_anas_after
    return
end

for i=1:N_tip
    if i==id, continue; end    
%     trace_i = reshape(traces(i,:,:),[frames 2]);
%     if min(vecnorm((trace_i-tips(id,:))'))<atm_min_dist
%         stop = true;
%         return
%     end
    % use the grid points 'occupied' by the vessels to test
    ix = get_vgp(tips(id,:));
    if vgp(i,ix(1),ix(2))
        stop = true;
        return
    end
end
end

function sc = check_self_intersection(id,v)
global tips vgp is_in_device %traces frames atm_min_dist

tip_new = tips(id)+v;
% trace_i = reshape(traces(id,:,:),[frames 2]);
% sc = min(vecnorm((trace_i-tip_new)'))<atm_min_dist;

if is_in_device(tip_new), sc = false; return; end
ix = get_vgp(tip_new);
sc = vgp(id,ix(1),ix(2));
end

function b = check_new_branch_at_tip(id)
% check if a new branch should be created at tip cell 'id'
global tips branch_tip_min branch_tip_max branch_tip_vegf

if branch_tip_max<=0, b = false; return; end

% if the tip cell has been stopped, do not branch
if check_tc_stop(id), b = false; return; end

vegf = vegf_conc(tips(id,1),tips(id,2));
prob = branch_tip_min + (branch_tip_max-branch_tip_min) ...
    * branch_tip_vegf*vegf / (1+branch_tip_vegf*vegf);

b = rand<prob;
end

function vd = get_vessel_density
global vgp sp %frames N_tip traces dx

% h = sp(2)-sp(1);
% vd = zeros(length(sp));
% for t=1:frames
%     for k=1:N_tip
%         if isnan(traces(k,t,1)), continue; end
%         ix = floor((traces(k,t,:)-sp(1))/h);
%         if t>1 && norm(reshape(traces(k,t,:)-traces(k,t-1,:),[1 2]))<dx
%             % if vessel has stopped movement, skip
%             continue
%         end
%         vd(ix(1),ix(2)) = vd(ix(1),ix(2)) + 1;
%     end
% end
% vd = vd/max(max(vd));

% use grid points 'occupied' by vessels. Much faster!
vd = reshape(sum(vgp,1),[length(sp) length(sp)]);
% vd = vd/max(max(vd));
end

function maintain_vgp
global N_tip tips vgp

for i=1:N_tip
    ix = get_vgp(tips(i,:));
    vgp(i,ix(1),ix(2)) = 1;
end
end

%% Tumor-related
function seed_tumor
global tumor_ctr N_tumor tumor_min_dist domain

tumor_ctr = repmat(domain*2,[N_tumor 1]); % initial out-of-domain positions
for i=1:N_tumor
    k = 1;
    while k<1e4
        tp_rejected = false;
        tp = gen_seed_pos;
        for j=1:N_tumor
            if i==j, continue; end
            if norm(tumor_ctr(j,:)-tp)<tumor_min_dist
                tp_rejected = true;
                break
            end
        end
        if ~tp_rejected, break; end
    end
    if tp_rejected, disp('seeding tumor failed'); return; end
    tumor_ctr(i,:) = tp;
end
end

function init_cont_field
global tumor_ctr N_tumor tumor_r sp tm nu nc ve

% tumor field
tm = zeros(length(sp));
for i=1:size(tm,1)
    x = sp(i);
    for j=1:size(tm,2)
        y = sp(j);
        for k=1:N_tumor
            r = (x-tumor_ctr(k,1))^2 + (y-tumor_ctr(k,2))^2;
            tm(i,j) = tm(i,j) + (1-tanh(r/tumor_r))/2;
        end
    end
end

% vegf field
ve = zeros(length(sp));
for i=1:size(ve,1)
    x = sp(i);
    for j=1:size(ve,2)
        y = sp(j);
        ve(i,j) = vegf_conc_init(x,y);
    end
end
ve = ve/max(max(ve)+1e-10);

% nutrient field
nu = nc*(1-tm/max(max(tm)+1e-10));
end

function f = field_only_in_device(f)
global sp is_in_device

% TODO: optimize
for i=1:length(sp)
    for j=1:length(sp)
        if ~is_in_device([sp(i) sp(j)]*0.98), f(i,j) = NaN; end
    end
end

% is_in_device = @(x) abs(x(1))+abs(x(2)/device_k) <= 1;

end

function update_cont_field
global sp tm nu ve tumor_D tumor_lam tumor_d tumor_dt ...
    nu_D nu_d nu_p nu_u nu_pd nc nc_v nu_dt ...
    vegf_D vegf_u vegf_d vegf_p0 vegf_p vegf_v0 vegf_dt

% get vessel density
vd = get_vessel_density;
vd = vd/max(max(vd));

% vegf equation:
% v_{,t} = D*div(grad(v)) - (u*y+d)*v + (p0+p*phi)*(v0-v)
% D = diffusivity, u = uptake by vessels (y), d = decay
% p0 = background production, p = production by tumors (phi)
% v0 = maximum vegf concentration
tmp = vegf_dt*vegf_D/(sp(2)-sp(1))^2;
src = (vegf_p0+vegf_p*tm).*(vegf_v0-ve) - (vegf_u*vd + vegf_d).*ve;
for i=2:size(ve)-1
    for j=2:size(ve)-1
        ve(i,j) = ve(i,j) + tmp*(ve(i+1,j)+ve(i-1,j)+ve(i,j+1)+ve(i,j-1)-4*ve(i,j)) ...
            + vegf_dt * src(i,j);
    end
end

% nutrient equation:
% n_{,t} = D*div(grad(n)) - (u*phi+d)*n + p*(1-phi)*(n0-n) + py*y*(n1-n)
% D = diffusivity, u = uptake, d = decay
% p = production by host (1-phi), n0 = nutrient concentration in host
% py = production by vessels (y), n1 = nutrient concentration in vessels
tmp = nu_dt*nu_D/(sp(2)-sp(1))^2;
src = nu_p*(1-tm).*(nc-nu) + nu_pd*vd.*(nc_v-nu) - (nu_u*tm + nu_d).*nu;
for i=2:size(nu)-1
    for j=2:size(nu)-1
        nu(i,j) = nu(i,j) + tmp*(nu(i+1,j)+nu(i-1,j)+nu(i,j+1)+nu(i,j-1)-4*nu(i,j)) ...
            + nu_dt * src(i,j);
    end
end

% tumor equation:
% phi_{,t} = M*div(grad(phi)) + (lambda*n - d)*phi
% M = mobility, lambda = mitosis rate, d = death rate
tmp = tumor_dt*tumor_D/(sp(2)-sp(1))^2;
src = (tumor_lam*nu - tumor_d) .* tm;
for i=2:size(tm)-1
    for j=2:size(tm)-1
        tm(i,j) = tm(i,j) + tmp*(tm(i+1,j)+tm(i-1,j)+tm(i,j+1)+tm(i,j-1)-4*tm(i,j)) ...
            + tumor_dt * src(i,j);
    end
end
end