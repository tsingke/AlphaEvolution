function [gbestx,gbestfitness,gbesthistory]=AE(N,D,ub,lb,MaxFEs,Func,FuncId)
% Algorithm Name: Alpha Evolution Algorithm (AE).
% gbestx: The global best solution ( gbestx = [x1,x2,...,xD]).
% gbestfitness: Record the fitness value of global best individual.
% gbesthistory: Record the history of changes in the global optimal fitness.
%---------------------------------------------------------------------------
%Initialization
FEs=0;

X=zeros(N,D);
R=X;
W=X;
newE=X;

f=inf(1,N);
gbesthistory=inf(1,MaxFEs);

%Building the candidate matrix
X=lb+(ub-lb).*rand(N,D);

for i=1:N
    f(i)=Func(X(i,:)',FuncId);
    FEs=FEs+1;
end

[gbestfitness,mi]=min(f);
gbestx=X(mi,:);
gbesthistory(1:FEs)=gbestfitness;

Pa=lb+(ub-lb).*rand(1,D);
Pb=lb+(ub-lb).*rand(1,D);
%---------------------------------------------------------------------------
while FEs<=MaxFEs
    %Sampling evolution matrix
    k=randi(N,N,1);
    E=X(k,:);
    [~,ind]=sort(f);
    R1=rand(N,D);R2=rand(N,D);S=randi([0,1],N,D);
    r=(ub-lb).*(2*R1.*R2-R2).*S;
    %Computing alpha
    alpha=exp(log(1-FEs/MaxFEs)-(4*(FEs/MaxFEs))^2);
    ar=alpha.*r;
    for i=1:N
        cab=FEs/MaxFEs;
        %Sampling W and L
        R(i,:)=X(ind(randi([1,length(1:find(k(i)==ind))])),:);
        W(i,:)=X(ind(randi([length(1:find(k(i)==ind)),N])),:);
        %Building evolutionary paths
        if rand<0.5
            A=X(randi(N,D,1),:);
            Pa=(1-cab)*Pa+cab*diag(A)';
            Ov=Pa;
        else
            K=ceil(N*rand);
            I1=[];I1=randperm(N,K);
            w=[];w=f(I1)./sum(f(I1));
            B=[];B=X(I1,:);
            Pb=(1-cab)*Pb+cab*(w*B);
            Ov=Pb;
        end
        I2=round(rand);
        sita=[];sita=I2*rand(1,D)+(1-I2)*rand*2;
        
        %AE Core operator
        newE(i,:)=Ov+ar(i,:)+sita.*(R(i,:)+E(i,:)-Ov-W(i,:));
        
        %Boundary constraint
        flagub=newE(i,:)>ub;flaglb=newE(i,:)<lb;
        newE(i,flagub)=(E(i,flagub)+ub)/2;newE(i,flaglb)=(E(i,flaglb)+lb)/2;
        newf=Func(newE(i,:)',FuncId);
        FEs=FEs+1;
        %Selection
        if newf<=f(k(i))
            f(k(i))=newf;
            X(k(i),:)=newE(i,:); 
            if f(k(i))<=gbestfitness
                gbestfitness=f(k(i));
                gbestx=X(k(i),:);
            end
        end
        gbesthistory(FEs)=gbestfitness;
    end
    if mod(FEs, MaxFEs/D)==0
        fprintf("AE, FEs: %d, fitess error = %e\n",FEs,gbestfitness);
    end
end % end while
%---------------------------------------------------------------------------
gbesthistory=gbesthistory(1:MaxFEs);
end