function GridModel
Par{1}=[1 3 30 80 0.9]; %Cost

N=50;
Par{2}=[N 1 0.5];%size of grid; neighbourhood  radius; nhood fraction for quorum
R = 1-ceil(sprand(N,N,0.001));
S = 1-ceil(sprand(N,N,0.001));%ones(N,N);%ceil(sprand(N,N,0.05))+1;
Par{3}=R;
Par{4}=S;
Par{5}=[400 10];%Generation runtime;  mutations per generation
Par{6}=[0 1 2]; %allele spectrum
Par{7}=[10 0];%rate of reading per competition; rate of image taking per generation - 0 means don't take
for i=1:20
Out{i}=RunGridModel(Par);
disp(i);

end
save('Run1','Par','Out')
return


function out=RunGridModel(Par)
%profile on

[Mr,Ms,Mc,M0,r]=getelements(Par{1});
[N,nh,fq]=getelements(Par{2});
nhood=2*nh+1;
nq=floor(fq*nhood.^2);

R=Par{3};
S=Par{4};
[Ng,Nmute]=getelements(Par{5});
AlleleSpec=Par{6};
Nas=length(AlleleSpec);

[Nvec,NIm]=getelements(Par{7});

I=ones(N,N);
I(1:N^2)=1:N^2;
I1=repmat(I,[3 3]);
Tor=I1((N-nh+(1:N+2*nh)),(N-nh+(1:N+2*nh))); %toroid conditions
Sf=S(Tor);
QS{1}=conv2(+(Sf==1),ones(nhood,nhood),'valid');
QS{2}=conv2(+(Sf==2),ones(nhood,nhood),'valid');
QR=(QS{1}>=nq).*(R==1)+(QS{2}>=nq).*(R==2);
QRf=QR(Tor);
PG=conv2(QRf,ones(nhood,nhood),'valid');

Cost=M0+Mr*(R>0)+Ms*(S>0)+Mc*QR;
Benefit=(PG>nq);%max(PG-nq,0)/nhood.^2;
Tcost=Cost.*(1-r*Benefit);
f=M0./Tcost;

%main loop
Nt=N.^2*Ng;


Z=3*(R)+S+1;
if NIm>0
Zv=zeros(N,N,round(Nt/Nim/N^2)+10);
Zv(:,:,1)=Z;
end
Nv=zeros(9,ceil(Nt/Nvec));
Nv(:,1)=hist(Z(:),1:9);
Nvs=Nv(:,1);
cn=zeros(1,Nt);

ns=1;
ic=0;
comp=1;
nvec=1;
tic
ndisp=1;
nmute=0;
while (ic<=Nt)
      [iw,jw,il,jl,cont]=competition(R,S,f);
      ic=ic+cont;
      comp=comp+1;
     
     
     Rold=R(il,jl);
     Sold=S(il,jl);
     
      Zold=3*(Rold)+Sold+1;
      Znew=3*(R(iw,jw))+S(iw,jw)+1;
      
      Nvs(Zold)=Nvs(Zold)-1;
      Nvs(Znew)=Nvs(Znew)+1;
    if min(Nvs)<0
    end
    
     R(il,jl)=R(iw,jw);
     S(il,jl)=S(iw,jw);
     
     %% Mutate
     
     Pmute=Nmute/N^2;
     if (rand(1)<cont*Pmute)
         nmute=nmute+1;
         i=ceil(rand(1)*N);
         j=ceil(rand(1)*N);
         al=ceil(rand(1)*2);
         ap=ceil(rand(1)*2);
         Zold=3*(R(i,j))+S(i,j)+1;
         if al==1
             if R(i,j)>0
                 MS=setdiff(AlleleSpec,R(i,j));
                 if Nas>2
                     MS=MS(ceil(rand(1)*length(MS)));
                 end
                 R(i,j)=MS;
             end
         else
             if S(i,j)>0
                 MS=setdiff(AlleleSpec,S(i,j));
                 if Nas>2
                     MS=MS(ceil(rand(1)*length(MS)));
                 end
                 S(i,j)=MS;
             end
         end
   
      Znew=3*(R(i,j))+S(i,j)+1;
      %Nv(:,comp)=Nv(:,comp-1);
      Nvs(Zold)=Nvs(Zold)-1;
      Nvs(Znew)=Nvs(Znew)+1;
     end
     
     if mod(comp,Nvec)==0
         nvec=nvec+1;
         cn(nvec)=ic;
         Nv(:,nvec)=Nvs;
     end
    %% Update
     I=mod(il+(-nh:nh)-1,N)+1;
     J=mod(jl+(-nh:nh)-1,N)+1;
     if Sold~=S(il,jl)
         if Sold>0
               QS{Sold}(I,J)=QS{Sold}(I,J)-1;
         end
         if S(il,jl)>0
             QS{S(il,jl)}(I,J)=QS{S(il,jl)}(I,J)+1;
         end
     end
     QR(I,J)=(QS{1}(I,J)>=nq).*(R(I,J)==1)+(QS{2}(I,J)>=nq).*(R(I,J)==2);
     %QRf=QR(Tor);
     I2=mod(il+(-2*nh:2*nh)-1,N)+1;
     J2=mod(jl+(-2*nh:2*nh)-1,N)+1;
     I3=mod(il+(-3*nh:3*nh)-1,N)+1;
     J3=mod(jl+(-3*nh:3*nh)-1,N)+1;

    PGs=conv2(QR(I3,J3),ones(nhood,nhood),'valid');
   Costs=M0+Mr*(R(I2,J2)>0)+Ms*(S(I2,J2)>0)+Mc*QR(I2,J2);
Benefits=PGs>nq;%max(PG(I,J)-nq,0)/nhood.^2;
Tcosts=Costs.*(1-r*Benefits);
f(I2,J2)=M0./Tcosts;
if NIm>0
    if floor(ic/Nim/N^2)>ns
        Zv(:,:,ns)=3*(R)+S+1;
        ns=ns+1;  
        disp(ic)
    end
end
if ic>N^2*ndisp
    disp(ic)
    ndisp=ndisp+1;
end
end
t=toc;
rm=find(cn>0,1,'last');
cn=cn(1:rm);
Nv=Nv(:,1:rm);
plot(cn,sum(Nv(1:3,:)),cn,Nv(5,:),cn,Nv(8,:),cn,Nv(9,:))
out{1}=cn;
out{2}=Nv;
if NIm>0
    out{3}=Zv;
end
return

function [iw,jw,il,jl,cont]=competition(R,S,f)
cont=0;
N=size(R,1);
dx=[-1 -1 0 1 1 1 0 -1];
dy=[0 1 1 1 0 -1 -1 -1];
while cont>-1
v=ceil(rand(1,2)*N);
i=v(1);
j=v(2);
in=ceil(rand(1)*8);
i1=mod(i+dx(in)-1,N)+1;
j1=mod(j+dy(in)-1,N)+1;
if R(i,j)+10*S(i,j)~=R(i1,j1)+10*S(i1,j1)
    if (f(i,j)+f(i1,j1))*rand(1)<f(i,j)
        iw=i;
        jw=j;
        il=i1;
        jl=j1;
    else
        iw=i1;
        jw=j1;
        il=i;
        jl=j;
    end   
    break
else
    cont=cont+1;
end
end
