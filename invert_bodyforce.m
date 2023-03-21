
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input Section

%load the mat file generated in build_bodyforce_Greens
%comment out line below if variables already loaded in workspace
load test_Greens



%specify range of beta values (weights on minimizing body forces)
%betas = linspace(6.1579,15,15);
betas = 10;

%number of realizations of strain rate for each beta value
num = 100;

%relative weight on fitting creep rate data (creeping faults)
Wcreep = 1;

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
    mhat = GL\d0;
    if ~isempty(PatchEnds)
        dhat = [G Gcreep([ind;ind],:)]*mhat;
    else
        dhat = G*mhat;
    end

 %propogate errors to strain rate
    if ~isempty(PatchEnds)
        Gsharp = (GG'*GG+beta^2*L'*L+Wcreep^2*Lcreep'*Lcreep)\GG';
    else
        Gsharp = (GG'*GG+beta^2*L'*L)\GG';
    end
    Cov_bf = Gsharp*Gsharp';  %note, data covariane is excluded because Gw is weighted

    Cov_exx = GExx*Cov_bf*GExx';
    Cov_exy = GExy*Cov_bf*GExy';
    Cov_eyy = GEyy*Cov_bf*GEyy';
     
     Exxs = GExx*mhat;
     Exys = GExy*mhat;
     Eyys = GEyy*mhat;
    
    
     Cov_Ve = GVe*Cov_bf*GVe';
     Cov_Vn = GVn*Cov_bf*GVn';
     
     Ves = GVe*mhat;
     Vns = GVn*mhat;
     
    Mhats(:,k) = mhat;
    
    
    %chi-squared
    chi_2(k) = sum((d./sig - dhat./sig).^2)/(length(d))
    
  
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

    
    
    save(savename,'chi_2','Mhats','Exx_realizations','Exy_realizations','Eyy_realizations','Ve_realizations','Vn_realizations')
    
    
    disp(['Completed ' num2str(k) ' of ' num2str(length(betas)) ' inversions'])
    
end

