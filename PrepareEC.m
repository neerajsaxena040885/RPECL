% -----------------------------------------------------------------
%      MIXED LOGIT MODEL FOR STATUS QUO ALTERNATIVE (SQ-RPL)
% -----------------------------------------------------------------

% Prepare a matrix of observed responses first
choice=zeros(NALTMAX,NCSMAX,NP);
for n=1:NP
    cs=XMAT(cp == n,2);
    yy=XMAT(cp == n,5);
    for t=1:NCSMAX
        choice(:,t,n)=yy((t-1)*NALTMAX+1:t*NALTMAX,1); 
    end
end
choice=repmat(choice,[1,1,1,NDRAWS]);

% Make draws for the Error Components SQ Alternative
disp('Creating Error Component draws for SQ Alternative.');
% Preparing the Master Draws matrix (MASTERDR) of dimensions 1 x NP x MASTERDRAWS
SQMASTERDR=ECSQmakedraws;   
SQMASTERDR=permute(SQMASTERDR,[3,2,1]);   %permute command is used to change the dimensions of a matrix... Read more on http://stackoverflow.com/questions/21100168/how-does-the-permute-function-in-matlab-work

% Define an empty smaller matrix to store draws from the Master Draws
SQDR=zeros(1,NP,NDRAWS);

% FIRST step is to extract a smaller matrix from Master Draws 
% Size of master Draws  -- NVxNPxMASTERDRAWS
% Size of smaller Draws -- NVxNPxNDRAWS
% LOGIC for extraction: Pick NDRAWS values at random from the Master Draws for every individual

for n=1:NP 
    % Select a value at random within the range [1,(AMPL-1)*NDRAWS]. This
    % marks the starting index for the NDRAWS to be selected. The Endig
    % index will be equal to the starting index + (NDRAWS-1). This way we
    % selected NDRAWS for an individual n
    if AMPL > 1
       startind=ceil(rand(1)*(MASTERDRAWS - NDRAWS - 1));
    else
       startind=1; 
    end
    endind=startind+(NDRAWS-1);
    SQDR(:,n,:)=SQMASTERDR(:,n,startind:endind);
end

% ECDR=repmat(SQDR,[1,1,1,NCSMAX]);
% ECDR01=reshape(ECDR,1,NCSMAX,NP,NDRAWS);
ECDR=reshape(SQDR,1,1,NP,NDRAWS);
ECDR01=repmat(ECDR,[1,NCSMAX,1,1]);

% clear SQDR ECDR SQMASTERDR


% Make draws for the Error Components Other Alternatives
disp('Creating Error Component draws for Other Alternatives.');
% Preparing the Master Draws matrix (MASTERDR) of dimensions 1 x NP x MASTERDRAWS
HYMASTERDR=ECHYmakedraws;   
HYMASTERDR=permute(HYMASTERDR,[4,3,2,1]);   %permute command is used to change the dimensions of a matrix... Read more on http://stackoverflow.com/questions/21100168/how-does-the-permute-function-in-matlab-work

% Define an empty smaller matrix to store draws from the Master Draws
HYDR=zeros((NALTMAX-1),NCSMAX,NP,NDRAWS);

% FIRST step is to extract a smaller matrix from Master Draws 
% Size of master Draws  -- NVxNPxMASTERDRAWS
% Size of smaller Draws -- NVxNPxNDRAWS
% LOGIC for extraction: Pick NDRAWS values at random from the Master Draws for every individual

for n=1:NP 
    % Select a value at random within the range [1,(AMPL-1)*NDRAWS]. This
    % marks the starting index for the NDRAWS to be selected. The Endig
    % index will be equal to the starting index + (NDRAWS-1). This way we
    % selected NDRAWS for an individual n
    if AMPL > 1
       startind=ceil(rand(1)*(MASTERDRAWS - NDRAWS - 1));
    else
       startind=1; 
    end
    endind=startind+(NDRAWS-1);
    HYDR(:,:,n,:)=HYMASTERDR(:,:,n,startind:endind);
end

% % ECDR=repmat(HYDR,[1,1,1,NCSMAX]);
% % ECDR02=reshape(ECDR,(NALTMAX-1),NCSMAX,NP,NDRAWS);
% ECDR=reshape(HYDR,(NALTMAX-1),1,NP,NDRAWS);
% ECDR02=repmat(ECDR,[1,NCSMAX,1,1]);
ECDR02=HYDR;

% clear HYDR ECDR HYMASTERDR

% Concatenate both ECDR01 and ECDR02 into a single matrix
ECDR=cat(1,ECDR01,ECDR02);  % Size of the Matrix will be NALTMAX X NCSMAX X NP X NDRAWS

clear ECDR01 ECDR02

% So now the Error Component matrix is ready.
% Matching it with the choice matrix and extracting the final matrix as (NALTMAX-1) X NCSMAX X NP X NDRAWS
% The resulting matrix, ECVECTOR, will be added to matrix v in the code "llgrad2.m"
ECVECTOR=zeros(NALTMAX-1,NCSMAX,NP,NMEM);

for draw=1:NDRAWS
    for n=1:NP
        for t=1:NCSMAX 
            if choice(1,t,n,draw)==1 && choice(2,t,n,draw)==0 && choice(3,t,n,draw)==0
               ECVECTOR(1,t,n,draw)=ECDR(2,t,n,draw)-ECDR(1,t,n,draw);
               ECVECTOR(2,t,n,draw)=ECDR(3,t,n,draw)-ECDR(1,t,n,draw);
            end
            if choice(1,t,n,draw)==0 && choice(2,t,n,draw)==1 && choice(3,t,n,draw)==0
               ECVECTOR(1,t,n,draw)=ECDR(1,t,n,draw)-ECDR(2,t,n,draw);
               ECVECTOR(2,t,n,draw)=ECDR(3,t,n,draw)-ECDR(2,t,n,draw);
            end
            if choice(1,t,n,draw)==0 && choice(2,t,n,draw)==0 && choice(3,t,n,draw)==1
               ECVECTOR(1,t,n,draw)=ECDR(1,t,n,draw)-ECDR(3,t,n,draw);
               ECVECTOR(2,t,n,draw)=ECDR(2,t,n,draw)-ECDR(3,t,n,draw);
            end
        end
    end
end

clear ECDR choice
