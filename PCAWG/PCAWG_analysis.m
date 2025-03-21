BKP.Properties.VarNames{8}='sv_type';
BKP.Properties.VarNames{9}='bkp1_type';
BKP.Properties.VarNames{10}='bkp2_type';
Sample=regexprep(BKP.Sample,'(.*):(.*)','$1');
SVC_ID=str2double(regexprep(BKP.Sample,'(.*):([0-9]*).*','$2'));
SVC_subID=1-logical(strcmp(regexprep(BKP.Sample,'(.*):([0-9]*).(.*)','$3'),''));

Junctions=sortrows([dataset(Sample,SVC_ID,'VarNames',{'Sample','SVC_id'}),BKP(:,2:end)],{'Sample','SVC_id','chr1','pos1'});
[Samples,IA,IC]=unique(Junctions.Sample);
Junctions.Sample=IC;

Junctions.sv_type=regexprep(Junctions.sv_type,'(.*):(.*)','$1');
Junctions.bkp1_type=regexprep(Junctions.bkp1_type,'(.*):(.*)','$1');
Junctions.bkp2_type=regexprep(Junctions.bkp2_type,'(.*):(.*)','$1');

STR1=Junctions.str1;
Junctions.str1(STR1=='+')='-';
Junctions.str1(STR1=='-')='+';
STR2=Junctions.str2;
Junctions.str2(STR2=='+')='-';
Junctions.str2(STR2=='-')='+';
Junctions=[dataset((1:1:length(Junctions))','VarNames','Idx'),Junctions];
save PCAWG_bkps.mat BKP Junctions Samples

%%
chrs=unique(Junctions.chr1);
Dist1=[];
Dist2=[];
Dist3=[];

Idx0=[];
Idx1=[];
Idx2=[];

INS=[];GAPS=[];PAR=[];

for si=1:length(Samples)
    junctions=Junctions(Junctions.Sample==si,:);
    junctions_2=junctions; 
        junctions_2.chr1=junctions.chr2;junctions_2.pos1=junctions.pos2;junctions_2.str1=junctions.str2;
        junctions_2.chr2=junctions.chr1;junctions_2.pos2=junctions.pos1;junctions_2.str2=junctions.str1;
        junctions_2.bkp1_type=junctions.bkp2_type;junctions_2.bkp2_type=junctions.bkp1_type;
    bkp=sortrows([junctions;junctions_2],{'chr1','pos1'});
    
    for ci=1:length(chrs)
        bkp_plus=bkp.pos1(strcmp(bkp.chr1,chrs{ci}) & bkp.str1=='+',:);
        bkp_minus=bkp.pos1(strcmp(bkp.chr1,chrs{ci}) & bkp.str1=='-',:);
        Dist3=[Dist3;diff(bkp_plus);diff(bkp_minus)];
        mat_dist=repmat(bkp_plus,1,length(bkp_minus))-repmat(bkp_minus',length(bkp_plus),1);
        mat_dist_1=mat_dist;mat_dist_1(mat_dist_1>0)=-3e8;
        dist1=-max(mat_dist_1,[],2);
        Dist1=[Dist1;dist1(dist1<3e8)];
        mat_dist_1=mat_dist;mat_dist_1(mat_dist_1<0)=3e8;
        dist2=min(mat_dist_1,[],1)';
        Dist2=[Dist2;dist2(dist2<3e8)];
    end
    
    % Short insertions or small overlaps
    idx=find(strcmp(bkp.chr1(1:end-1),bkp.chr1(2:end)) & bkp.str1(1:end-1)=='+' & bkp.str1(2:end)=='-' & diff(bkp.pos1)<=20000); 
    while length(idx)>0
        INS=[INS;bkp(unique([idx;idx+1]),:)];
        Idx0=[Idx0;unique(bkp.Idx(unique([idx;idx+1])))];
        % excluding short insertions/small overlapping ends from adjacent
        % breakpoint tally
        bkp(unique([idx;idx+1]),:)=[];
        idx=find(strcmp(bkp.chr1(1:end-1),bkp.chr1(2:end)) & bkp.str1(1:end-1)=='+' & bkp.str1(2:end)=='-' & diff(bkp.pos1)<=20000);
    end
   
    % Adjacent gapped breakpoints
    idx=find(strcmp(bkp.chr1(1:end-1),bkp.chr1(2:end)) & bkp.str1(1:end-1)=='-' & bkp.str1(2:end)=='+' & diff(bkp.pos1)<=20000);
    GAPS=[GAPS;bkp(unique([idx;idx+1]),:)];
    Idx1=[Idx1;unique(bkp.Idx(unique([idx;idx+1])))];
   
    % Adjacent parallel breakpoints
    idx=find(strcmp(bkp.chr1(1:end-1),bkp.chr1(2:end)) & bkp.str1(1:end-1)==bkp.str1(2:end) & diff(bkp.pos1)<=20000);
    PAR=[PAR;bkp(unique([idx;idx+1]),:)];
    Idx2=[Idx2;unique(bkp.Idx(unique([idx;idx+1])))];
end

% These are single breakpoints
Idx3=setdiff(1:1:length(Junctions),unique([Idx0;Idx1;Idx2]));
SingleJunctions=Junctions(Idx3,:);

%% Tally of different types of breakpoints
bkp_tally=[];
for si=1:length(Samples)
    junctions=Junctions(Junctions.Sample==si,:);
    sj=SingleJunctions(SingleJunctions.Sample==si,:);
    ins=INS(INS.Sample==si,:);par=PAR(PAR.Sample==si,:);gaps=GAPS(GAPS.Sample==si,:);
    svs=unique(junctions.SVC_id);
    for ci=1:length(svs)
        bkp_tally=[bkp_tally;si,ci,sum(junctions.SVC_id==svs(ci)),sum(sj.SVC_id==svs(ci)),sum(ins.Sample==si & ins.SVC_id==svs(ci)),sum(gaps.Sample==si & gaps.SVC_id==svs(ci)),sum(par.Sample==si & par.SVC_id==svs(ci))];
    end
end
SV_tally=dataset(bkp_tally(:,1),bkp_tally(:,2),bkp_tally(:,3),bkp_tally(:,4),bkp_tally(:,5),bkp_tally(:,6),bkp_tally(:,7),...
    'VarNames',{'Sample','SVC_id','total','single_bkps','ins_bkps','gap_bkps','par_bkps'});

%% Groups of adjacent parallel breakpoints
pos_diff=[diff(PAR.pos1);-1];
pos_diff(diff(PAR.Sample)>0 | ~strcmp(PAR.chr1(1:end-1),PAR.chr1(2:end)) | PAR.str1(1:end-1)~=PAR.str1(2:end))=-1;
pos_diff(pos_diff>2e4)=-1;

s_idx=find(pos_diff==-1);
adj_idx=[];
adj_idx(1:s_idx(1))=1;
for i=1:length(s_idx)-1
    adj_idx(s_idx(i)+1:s_idx(i+1))=i+1;
end

PAR=[dataset(adj_idx','VarNames',{'par_id'}),PAR];

%%
save PCAWG_analysis_results.mat Junctions Idx0 Idx1 Idx2 Idx3 Dist1 Dist2 Dist3 INS GAPS PAR SV_tally

%% Adjacent parallel breakpoints with adjacent cis breakpoints

[Shared,IA,IB]=intersect(PAR(:,2:7),GAPS(:,1:6));
i=1;
A=[];B=[];
aid=1;bid=1;
while i<max(PAR.par_id)
    curr_par=PAR(PAR.par_id==i,:);
    if curr_par.str1(1)=='+'
        curr_bkp=curr_par(1,:);
        [~,~,cid]=intersect(curr_bkp(1,2:7),Shared);
        if ~isempty(cid)
            add=[GAPS(IB(cid)-1,:);curr_par(:,2:end)];
            A=[A;dataset(repmat(aid,size(add(:,1))),'VarNames',{'cid'}),add];
            aid=aid+1;
        end
        i=i+1;
    else
        curr_bkp=curr_par(end,:);
        i=i+1;
        next_par=PAR(PAR.par_id==i,:);
        next_bkp=next_par(1,:);
        if curr_bkp.Sample==next_bkp.Sample & strcmp(curr_bkp.chr1,next_bkp.chr1) & next_bkp.str1=='+' & next_bkp.pos1-curr_bkp.pos1<=2e4
            B=[B;dataset(repmat(bid,size(curr_par(:,1))),'VarNames',{'cid'}),curr_par(:,2:end);dataset(repmat(bid,size(next_par(:,1))),'VarNames',{'cid'}),next_par(:,2:end)];
            bid=bid+1;
            i=i+1;
        else
            [~,~,cid]=intersect(curr_bkp(1,2:7),Shared);
            if ~isempty(cid)
                add=[curr_par(:,2:end);GAPS(IB(cid)+1,:)];
                A=[A;dataset(repmat(aid,size(add(:,1))),'VarNames',{'cid'}),add];
                aid=aid+1;
            end
        end
    end
end

PAR_with_GAP=A;
PAR_with_PAR=B;

save adj_with_cis_bkps.mat PAR*

%% Landscape of breakpoints on each chromosome arm
load hg19.mat
clear A B C D E F G H
for i=1:length(ChrArms)
    currArm=ChrArms(i,:);
    currChr=ChrArmsLabels{i}(1:end-1);
    Junctioni=Junctions((strcmp(Junctions.chr1,currChr) & Junctions.pos1>=currArm.left & Junctions.pos1<=currArm.right)...
                      | (strcmp(Junctions.chr2,currChr) & Junctions.pos2>=currArm.left & Junctions.pos2<=currArm.right),:);
                  
    SJ=SingleJunctions((strcmp(SingleJunctions.chr1,currChr) & SingleJunctions.pos1>=currArm.left & SingleJunctions.pos1<=currArm.right)...
                      |(strcmp(SingleJunctions.chr2,currChr) & SingleJunctions.pos2>=currArm.left & SingleJunctions.pos2<=currArm.right),:);
    
    PJ=PAR(strcmp(PAR.chr1,currChr) & PAR.pos1>=currArm.left & PAR.pos1<=currArm.right,:);
    INSJ=INS(strcmp(INS.chr1,currChr) & INS.pos1>=currArm.left & INS.pos1<=currArm.right,:);
    GAPJ=GAPS(strcmp(GAPS.chr1,currChr) & GAPS.pos1>=currArm.left & GAPS.pos1<=currArm.right,:);
    
    for si=1:length(Samples)
        A(si,i)=sum(SJ.Sample==si); %Junctions with breakpoints that do not have adjacent breakpoints
        B(si,i)=sum(PJ.Sample==si); %Junctions with breakpoints that have an adjacent parallel breakpoint
        C(si,i)=sum(INSJ.Sample==si); % Junctions with breakpoints that are insertion/overlapping adjacent breakpoints
        D(si,i)=sum(GAPJ.Sample==si); % Junctions with breakpoints that have an adjacent gapped breakpoint
        junctions=Junctioni(Junctioni.Sample==si,:);
        junctions_2=junctions; 
        junctions_2.chr1=junctions.chr2;junctions_2.pos1=junctions.pos2;junctions_2.str1=junctions.str2;
        junctions_2.chr2=junctions.chr1;junctions_2.pos2=junctions.pos1;junctions_2.str2=junctions.str1;
        junctions_2.bkp1_type=junctions.bkp2_type;junctions_2.bkp2_type=junctions.bkp1_type;
        bkp=sortrows([junctions;junctions_2],{'chr1','pos1'});
        bkp=bkp(strcmp(bkp.chr1,currChr) & bkp.pos1>=currArm.left & bkp.pos1<=currArm.right,:);
        all_adj=sortrows(unique([PJ(PJ.Sample==si,2:10);INSJ(INSJ.Sample==si,1:9);GAPJ(GAPJ.Sample==si,1:9)]),{'chr1','pos1'});
        E(si,i)=length(bkp(:,1:9)); % All breakpoints
        F(si,i)=length(setdiff(bkp(:,1:9),all_adj)); % All breakpoints that do not have adjacent breakpoints; their partners may have adjacent breakpoints
        G(si,i)=length(intersect(PJ(PJ.Sample==si,2:10),GAPJ(GAPJ.Sample==si,1:9)));
    end
end
bkp_density=E./repmat((ChrArms.right-ChrArms.left)',size(E(:,1)));
idx=find(bkp_density>=1e-6 & E>=100);
bkpd=bkp_density(idx);
I=sortrows([A(idx),B(idx),C(idx),E(idx),bkpd],4);

figure
subplot(2,1,1),
bar([I(:,1:3),I(:,4)-I(:,1)-I(:,2)-I(:,3)],'stacked')
subplot(2,1,2),
bar(I(:,5)*1e6,1)
%% chr4p breakpoint distribution
bins=0:1e5:5e7;
% Single breakpoints from all samples
pos=double(sort([SingleJunctions.pos1(strcmp(SingleJunctions.chr1,'4'));SingleJunctions.pos2(strcmp(SingleJunctions.chr2,'4'))]));
counts=histc(pos,bins);
figure,subplot(2,1,1)
plot((bins(1:end-1)+bins(2:end))/2/1e6,counts(1:end-1)/mean(counts(1:end-1)),'o')
box off
set(gca,'TickDir','out','FontSize',16)
subplot(2,1,2)
plot(pos(1:end-1)/1e6,diff(pos)/mean(diff(pos)),'o'),box off
set(gca,'Ylim',[0.01 100],'Xlim',[0 50])
set(gca,'YScale','log','TickDir','out','FontSize',16)

figure,
diff_pos=diff(pos(pos<5e7));
histogram(log10(diff(pos)/mean(diff(pos))),-5:0.15:2)

%Breakpoints from all samples
pos=double(sort([Junctions.pos1(strcmp(Junctions.chr1,'4'));Junctions.pos2(strcmp(Junctions.chr2,'4'))]));
counts=histc(pos,bins);
figure,subplot(2,1,1)
plot((bins(1:end-1)+bins(2:end))/2/1e6,counts(1:end-1)/mean(counts(1:end-1)),'o')
box off
set(gca,'TickDir','out','FontSize',16)
subplot(2,1,2)
plot(pos(1:end-1)/1e6,diff(pos)/mean(diff(pos)),'o'),box off
set(gca,'Ylim',[0.01 100],'Xlim',[0 50])
set(gca,'YScale','log','TickDir','out','FontSize',16)

%Adjacent Parallel Breakpoints
pos=double(sort(PAR.pos1(strcmp(PAR.chr1,'4'))));
counts=histc(pos,bins);
figure,subplot(2,1,1)
plot((bins(1:end-1)+bins(2:end))/2/1e6,counts(1:end-1)/mean(counts(1:end-1)),'o')
box off
set(gca,'TickDir','out','FontSize',16)
subplot(2,1,2)
plot(pos(1:end-1)/1e6,diff(pos)/mean(diff(pos)),'o'),box off
set(gca,'Ylim',[0.01 100],'Xlim',[0 50])
set(gca,'YScale','log','TickDir','out','FontSize',16)

%Insertions
pos=double(sort(INS.pos1(strcmp(PAR.chr1,'4'))));
counts=histc(pos,bins);
figure,subplot(2,1,1)
plot((bins(1:end-1)+bins(2:end))/2/1e6,counts(1:end-1)/mean(counts(1:end-1)),'o')
box off
set(gca,'TickDir','out','FontSize',16)
subplot(2,1,2)
plot(pos(1:end-1)/1e6,diff(pos)/mean(diff(pos)),'o'),box off
set(gca,'Ylim',[0.01 100],'Xlim',[0 50])
set(gca,'YScale','log','TickDir','out','FontSize',16)

