
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input Section

%specify range of beta values (weights on minimizing body forces)
betas = [40 45 50];

%compute uncertainties?  set ucertainty = true or false
%computing uncertainties is slow
uncertainty = false;
%number of realizations of strain rate for each beta value
num = 500;

%relative weight on fitting creep rate data (creeping faults)
Wcreep = 1;

% optional two-step minimization of strain rates below threshold value
% set twostep = true or false
% relative weight (gamma) on minimizing strain rates below strain_threshold
%(micro-strain per year)
twostep = false;
gamma = 400;
strain_threshold = 9e-3;  %micro-strain per year
%x and y boundaries for strain rate minimization (strain rates not
%minimized outside of these boundaries)
minimize_xbound = [-600 1000];
minimize_ybound = [-800 1000];

%name of mat file for saving inversion results
savename = 'test_inversion';

%% No need to modify below
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% find cells with data
ind = ~isnan(Ve_centroids);

G = [Ge_x(ind,:) Ge_y(ind,:);Gn_x(ind,:) Gn_y(ind,:)];
d = [Ve_centroids(ind);Vn_centroids(ind)];

Sige = Sige_centroids(ind);
Sign = Sign_centroids(ind);

%build creep matrix
if ~isempty(PatchEnds)
    build_creepG
end

%% regularization weights

%build_triangle_weights -- L, to penalize magnitude of body forces
I = eye(size(Ge_x,2));
%I = diag(w_tri);
%I = con_tri.*repmat(w_tri,1,size(con_tri,2));
L = [I 0*I;0*I I];

if ~isempty(PatchEnds)
    L = [L zeros(size(L,1),length(PatchCreepRates))];
    Lcreep = [zeros(length(PatchCreepRates),size(G,2)) eye(length(PatchCreepRates))];
end




sig = [Sige;Sign];


   
%matrix for computing strain rates
if ~isempty(PatchEnds)
    GExx = [GExx_x GExx_y GExx_creep];
    GExy = [GExy_x GExy_y GExy_creep];
    GEyy = [GEyy_x GEyy_y GEyy_creep];
else
    GExx = [GExx_x GExx_y];
    GExy = [GExy_x GExy_y];
    GEyy = [GEyy_x GEyy_y];
end


%velocities
if ~isempty(PatchEnds)
    GVe = [Ge_x Ge_y Geast_creep];
    GVn = [Gn_x Gn_y Gnorth_creep];
else
    GVe = [Ge_x Ge_y];
    GVn = [Gn_x Gn_y];   
end

%build realizations
Exx_realizations = [];
Exy_realizations = [];
Eyy_realizations = [];

Ve_realizations = [];
Vn_realizations = [];

if ~isempty(PatchEnds)
    GG = [G Gcreep([ind;ind],:)];
else
    GG = G; 
end
GG = GG./repmat(sig,1,size(GG,2));

dd = d./sig;

for k=1:length(betas)

    disp(['Beginning ' num2str(k) ' of ' num2str(length(betas)) ' inversions'])
    
    beta = betas(k);

    if ~isempty(PatchEnds)
        GL = sparse([GG;[beta*L;Wcreep*Lcreep]]);
        d0 = sparse([dd;zeros(size(L,1),1);Wcreep*PatchCreepRates]);
    else
        GL = sparse([GG;beta*L]);
        d0 = sparse([dd;zeros(size(L,1),1)]);
    end

    tic; mhat0 = GL\d0; T = toc;
    
     
    disp(['Completed first inversion step in ' num2str(T) ' seconds.']);

     if twostep 
    
         Exxs = GExx*mhat0;
         Exys = GExy*mhat0;
         Eyys = GEyy*mhat0;
         max_shear = sqrt((Exxs-Eyys).^2 + Exys.^2);
    
         ind_boundaries = tri_centroids(:,1)>minimize_xbound(1) & tri_centroids(:,1)<minimize_xbound(2) & tri_centroids(:,2)>minimize_ybound(1) & tri_centroids(:,2)<minimize_ybound(2);
         ind_threshold = (max_shear < strain_threshold) & ind_boundaries;

         disp(['Minimizing strain rates in ' num2str(sum(ind_threshold)/length(ind_threshold)*100) ' % of cells.'])

         L_iter = [GExx(ind_threshold,:);GExy(ind_threshold,:);GEyy(ind_threshold,:)];
         GL_iterate = [GL; gamma*L_iter] ;
         d0_iterate = [d0; zeros(3*sum(ind_threshold),1)];


         tic;  mhat = lsqr(GL_iterate,d0_iterate,[],1500,[],[],mhat0);  T = toc;
         disp(['Completed second (minimization) inversion step in ' num2str(T) ' seconds.']);


     else

         mhat = mhat0;  %use first step as solution

     end

    %predicted velocities using mhat computed above
    if ~isempty(PatchEnds)
        dhat = [G Gcreep([ind;ind],:)]*mhat;
    else
        dhat = G*mhat;
    end


    %propogate errors to strain rate
    if uncertainty 

        if twostep
    
           if ~isempty(PatchEnds)
               Ginv = (GG'*GG+beta^2*L'*L + Wcreep^2*Lcreep'*Lcreep + gamma*L_iter'*L_iter);
           else
               Ginv = (GG'*GG+beta^2*L'*L + gamma*L_iter'*L_iter);
           end
    
        else
    
           if ~isempty(PatchEnds)
               Ginv = (GG'*GG+beta^2*L'*L + Wcreep^2*Lcreep'*Lcreep);
           else
               Ginv = (GG'*GG+beta^2*L'*L);
           end
    
        end
    
        Gsharp = Ginv\GG';
       
        %error propagation for strain rates
        Cov_bf = Gsharp*Gsharp';  %note, data covariance is excluded because Gw is weighted
        Cov_exx = GExx*Cov_bf*GExx';
        Cov_exy = GExy*Cov_bf*GExy';
        Cov_eyy = GEyy*Cov_bf*GEyy';
    
        %error propagation for velocities
        Cov_Ve = GVe*Cov_bf*GVe';
        Cov_Vn = GVn*Cov_bf*GVn';

    end

    %computed strain rates and velocities
    Exxs = GExx*mhat;
    Exys = GExy*mhat;
    Eyys = GEyy*mhat;
     
    Ves = GVe*mhat;
    Vns = GVn*mhat;
     
    Mhats(:,k) = mhat;  %store all solutions
    
    
    %chi-squared
    chi_2(k) = sum((d./sig - dhat./sig).^2)/(length(d))
    

    if uncertainty
  
        %draw realizations from multivariate Gaussian distribution
        tic 
    
        J = size(Cov_exx,1);
          
        noise = mvnrnd(zeros(1,J),full(Cov_exx),num);
        Exx_realizations = [Exx_realizations repmat(Exxs,1,num) + noise'];
        
         noise = mvnrnd(zeros(1,J),full(Cov_exy),num);
         Exy_realizations = [Exy_realizations repmat(Exys,1,num) + noise'];
       
         noise = mvnrnd(zeros(1,J),full(Cov_eyy),num);
        Eyy_realizations = [Eyy_realizations repmat(Eyys,1,num) + noise'];
    
        
        noise = mvnrnd(zeros(1,J),full(Cov_Ve),num);
        Ve_realizations = [Ve_realizations repmat(Ves,1,num) + noise'];
        noise = mvnrnd(zeros(1,J),full(Cov_Vn),num);
        Vn_realizations = [Vn_realizations repmat(Vns,1,num) + noise'];
    
        T = toc;
    
        disp(['Completed ' num2str(num) ' realizations in ' num2str(T) ' seconds.']);

    else

        Exx_realizations = [Exx_realizations Exxs];
        Exy_realizations = [Exy_realizations Exys];
        Eyy_realizations = [Eyy_realizations Eyys];

        Ve_realizations = [Ve_realizations Ves];
        Vn_realizations = [Vn_realizations Vns];
    end

    
    save(savename,'chi_2','Mhats','Exx_realizations','Exy_realizations','Eyy_realizations','Ve_realizations','Vn_realizations')
    
    
    disp(['Completed ' num2str(k) ' of ' num2str(length(betas)) ' inversions'])
    
end

