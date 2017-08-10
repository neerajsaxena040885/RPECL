% -----------------------------------------------------------------
%       MIXED LOGIT MODEL FOR STATUS QUO ALTERNATIVE (SQ-RPL)
% -----------------------------------------------------------------

% Array of draws has dimension MASTERDRAWS x NP x 1.

% -----------------------------------------------------------------
%                            CHANGE LOG

% 12/06/2016: Introducing other drawing techniques as well. 
% 12/06/2016: Amplifying the size of draw matrix to [MASTERDRAW x NP x
% 1]. Call this as the Master Draws matrix 
% -----------------------------------------------------------------

function dr=ECSQmakedraws
 
global NV NP DRAWTYPE MASTERDRAWS NALTMAX

% Define an empty array to store the draws
dr=zeros(MASTERDRAWS,NP,1);


% ----- HALTON DRAWS (DRAWTYPE=2) & SHIFTED HALTON DRAWS (DRAWTYPE=3) -----
%
if DRAWTYPE == 2 || DRAWTYPE == 3   % Halton draws
   h=primes(100);          %Generates prime numbers below 100. The array is of size 1X25 in this case. There are 25 prime numbers below 100      % Must create for all people together
   k=1;
   while size(h,2) < (NV+NALTMAX)       %size(h,2) returns the number of prime numbers. Since we have to maintain separate prime for each random coeff, thus checking if the number of primes is sufficient, i.e. Number of primes > Number of parameters to be estimated
       h=primes(k.*100);
       k=k+1;
   end
   h=h(1,NV+1);
   for j=1:1
       hh=h(1,j);   % Selecting column j from prime number array as the prime for variable j
       draws=[0];   % The first element of draw array is set to zero initially. This element can be dropped later
       test=0;
       b=1;
       while test == 0
            drawsold=draws;
            for m=1:(hh-1);
                dd=m./(hh.^b);
                draws=[draws ; drawsold + dd];
                test=size(draws,1) >= ((NP.*MASTERDRAWS) + 10);    % +10 here means that we are generating extra 10 draws, which means I would be dropping first 10 elements in next steps. This is being done in order to minimise correlation 
                if test == 1
                   break
                end
            end
            b=b+1;    
       end
       draws=draws(11:(NP.*MASTERDRAWS)+10,1);   % Dropping the first 10 elements from further processing. This is being done in order to minimise correlation among prime numbers
       
       if DRAWTYPE == 3     % Shifted and Shuffled Halton Draws
            draws=draws+rand(1,1);     %rand(1,1) creates a cell of randomly generated number between 0 and 1      %Shift: one shift for entire sequence
            draws=draws-floor(draws);  %floor function to round down the value. Two cases can arise: either draws=draws OR draws=draws-1... This shift reduces the level of correlation
            draws=reshape(draws,MASTERDRAWS,NP);
            for n=1:NP                  %Shuffle for each person separately
               rr=rand(MASTERDRAWS,1);       % generate a column array of uniformly distributed numbers between (0,1)of size MASTERDRAWS
               [rr rrid]=sort(rr);      % sort the rr array and store the values in rr array and index of sorted elements in rrid
               draws(:,n)=draws(rrid,n);    % Shuffle the draws column for nth person based on the rrid generated above. This considerably brings down the level of correlation
            end;
            draws=reshape(draws,NP.*MASTERDRAWS,1);  % a column array draws having NP*MASTERDRAWS rows
       end
       % Since Error Components are Standard Normally distributed.
       draws=-sqrt(2).*erfcinv(2.*draws);  %Take inverse cum normal
       dr(:,:,j)=reshape(draws,MASTERDRAWS,NP);  % a 3-D matrix of shifted and shuffled Halton draws
    end
end


% ------------------ MLHS DRAWS (DRAWTYPE=4) -------------------
%
if DRAWTYPE == 4
    h=0:(MASTERDRAWS-1);
    h=h'./MASTERDRAWS;
    for j=1:1
        for n=1:NP
            draws=h+rand(1,1)./MASTERDRAWS;    % rand(1,1) is divided by MASTERDRAWS and the answer is added to h      %Shift: Different shift for each person
            rr=rand(MASTERDRAWS,1);
            [rr rrid]=sort(rr);
            draws=draws(rrid,1);          %Shuffle
            % Since Error Components are Standard Normally distributed.
            draws=-sqrt(2).*erfcinv(2.*draws);  %Take inverse cum normal
            dr(:,n,j)=draws;
        end
    end
end
  