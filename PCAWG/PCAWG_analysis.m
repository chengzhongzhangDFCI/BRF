temp=textscan(fopen('all_sv_classifications.trimmed.txt','r'),'%s\t%d\t%d\t%s\t%d\t%d\t%c\t%c\t%s\t%s\t%s\t%s\t\n');
BKP=dataset(temp{9},temp{1},temp{2},temp{7},temp{4},temp{5},temp{8},temp{10},temp{11},temp{12},'VarNames',...
           {'Sample','chr1','pos1','str1','chr2','pos2','str2','sv_type','bkp1_type','bkp2_type'});
Sample=regexprep(BKP.Sample,'(.*):(.*)','$1');
SVC_ID=str2double(regexprep(BKP.Sample,'(.*):([0-9]*).*','$2'));
%SVC_subID=1-logical(strcmp(regexprep(BKP.Sample,'(.*):([0-9]*).(.*)','$3'),''));

Junctions=sortrows([dataset(Sample,SVC_ID,'VarNames',{'Sample','SVC_id'}),BKP(:,2:end)],{'Sample','SVC_id','chr1','pos1'});
[Samples,IA,IC]=unique(Junctions.Sample);
Junctions.Sample=IC;

Junctions.sv_type=regexprep(Junctions.sv_type,'(.*):(.*)','$1'); % only consider type of SV events
Junctions.bkp1_type=regexprep(Junctions.bkp1_type,'(.*):(.*)','$1');
Junctions.bkp2_type=regexprep(Junctions.bkp2_type,'(.*):(.*)','$1');

STR1=Junctions.str1;
Junctions.str1(STR1=='+')='-'; % revert to the convention used in the paper
Junctions.str1(STR1=='-')='+';
STR2=Junctions.str2;
Junctions.str2(STR2=='+')='-';
Junctions.str2(STR2=='-')='+';
Junctions=[dataset((1:1:length(Junctions))','VarNames','Idx'),Junctions]; % Each row is a junction between two breakpoints 
%save PCAWG_bkps.mat BKP Junctions Samples

%% Find different classes of adjacent breakpoints
INS=[];GAPS=[];PAR=[]; % Insertion breakpoints, gapped breakpoints, parallel breakpoints that are adjacent
NonINSBKPS=[]; % Breakpoints that are not related to insertions

chrs=unique(Junctions.chr1);
dist_threshold=20000; % Threshold for adjacent breakpoints

%Dist1=[];
%Dist2=[];
%Dist3=[];

Idx0=[];
Idx1=[];
Idx2=[];

for si=1:length(Samples) % process for each individual sample
    junctions=Junctions(Junctions.Sample==si,:);
    junctions_2=junctions; 
        junctions_2.chr1=junctions.chr2;junctions_2.pos1=junctions.pos2;junctions_2.str1=junctions.str2; 
        junctions_2.chr2=junctions.chr1;junctions_2.pos2=junctions.pos1;junctions_2.str2=junctions.str1;
        junctions_2.bkp1_type=junctions.bkp2_type;junctions_2.bkp2_type=junctions.bkp1_type;
    % collect and sort breakpoints from all the junctions by coordinate
    bkp=sortrows([junctions;junctions_2],{'chr1','pos1'});
    
    % Generate breakpoint distributions
%     for ci=1:length(chrs)
%         bkp_plus=bkp.pos1(strcmp(bkp.chr1,chrs{ci}) & bkp.str1=='+',:); % (+) breakpoints
%         bkp_minus=bkp.pos1(strcmp(bkp.chr1,chrs{ci}) & bkp.str1=='-',:); % (-) breakpoints
%         Dist3=[Dist3;diff(bkp_plus);diff(bkp_minus)]; % distance between nearest breakpoints (+/+ or -/-)
%         mat_dist=repmat(bkp_plus,1,length(bkp_minus))-repmat(bkp_minus',length(bkp_plus),1); % pairwise distance between each (+) breakpoint and each (-) breakpoint 
%         mat_dist_1=mat_dist;mat_dist_1(mat_dist_1>0)=-3e8; 
%         dist1=-max(mat_dist_1,[],2); % distance between nearest (+) and (-) breakpoints (for insertions/overlapping breakpoints)
%         Dist1=[Dist1;dist1(dist1<3e8)];
%         mat_dist_1=mat_dist;mat_dist_1(mat_dist_1<0)=3e8; % distance between nearest (-) and (+) breakpoints (gapped)
%         dist2=min(mat_dist_1,[],1)';
%         Dist2=[Dist2;dist2(dist2<3e8)];
%     end
    
    % Short insertions or small overlaps (adjacent (+) and (-) breakpoints)
    idx=find(strcmp(bkp.chr1(1:end-1),bkp.chr1(2:end)) & bkp.str1(1:end-1)=='+' & bkp.str1(2:end)=='-' & diff(bkp.pos1)<=dist_threshold); 
    while length(idx)>0
        INS=[INS;bkp(unique([idx;idx+1]),:)]; % insertion breakpoints
        Idx0=[Idx0;unique(bkp.Idx(unique([idx;idx+1])))]; % indices of insertion breakpoints
        bkp(unique([idx;idx+1]),:)=[]; % iteratively remove short insertions/small overlapping ends from adjacent breakpoint tally
        idx=find(strcmp(bkp.chr1(1:end-1),bkp.chr1(2:end)) & bkp.str1(1:end-1)=='+' & bkp.str1(2:end)=='-' & diff(bkp.pos1)<=dist_threshold);
    end
    NonINSBKPS=[NonINSBKPS;bkp]; % breakpoints without an adjacent breakpoint that can possibly make an insertion
    
    % Adjacent gapped breakpoints
    idx=find(strcmp(bkp.chr1(1:end-1),bkp.chr1(2:end)) & bkp.str1(1:end-1)=='-' & bkp.str1(2:end)=='+' & diff(bkp.pos1)<=dist_threshold);
    GAPS=[GAPS;bkp(unique([idx;idx+1]),:)];
    Idx1=[Idx1;unique(bkp.Idx(unique([idx;idx+1])))];
   
    % Adjacent parallel breakpoints
    idx=find(strcmp(bkp.chr1(1:end-1),bkp.chr1(2:end)) & bkp.str1(1:end-1)==bkp.str1(2:end) & diff(bkp.pos1)<=dist_threshold);
    PAR=[PAR;bkp(unique([idx;idx+1]),:)];
    Idx2=[Idx2;unique(bkp.Idx(unique([idx;idx+1])))];
end

% These are single breakpoints: neither these nor their partners have
% adjacent breakpoints
Idx3=setdiff(1:1:length(Junctions),unique([Idx0;Idx1;Idx2]));
SingleJunctions=Junctions(Idx3,:);

%save NonINS.mat NonINSBKPS

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
max(PAR.par_id); % number of adjacent parallel breakpoints (can have trios etc.)
length(unique(PAR.Sample)); % number of samples with adjacent parallel breakpoints;
length(unique(PAR.Idx)); % number of junctions where one or both breakpoints have an adjacent parallel breakpoint
tmp=sortrows(PAR,{'Idx','par_id'});
sum(diff(tmp.Idx)==0 & diff(tmp.par_id)==0); % foldback junctions: both breakpoints (with the same Idx) belong to a single PAR group (par_id)

%%
%save PCAWG_analysis_results.mat Junctions Idx0 Idx1 Idx2 Idx3 Dist1 Dist2 Dist3 INS GAPS PAR SV_tally

%% Estimation of probability of independent breakpoint generation
PAR_loci=[];
load hg19.mat
for si=1:length(Samples)
    par=PAR(PAR.Sample==si,:);
    noninsbkps=NonINSBKPS(NonINSBKPS.Sample==si,:);
    par_ids=unique(par.par_id);
    for pid=1:length(par_ids)
        par_id=par_ids(pid);
        no_bkp=sum(par.par_id==par_id);
        bkp1=par(find(par.par_id==par_id,1,'first'),:); % range of adjacent parallel breakpoints
        bkp2=par(find(par.par_id==par_id,1,'last'),:);  % in each group
        chr=bkp1.chr1; % all breakpoints are on the same chromosome
        str=bkp1.str1; % and have the same orientation
        if str=='-'
            % for (-) breakpoints, find the (+) breakpoint to the left of
            % the first (-) breakpoint in the group
            adjbkp_pos=max(noninsbkps.pos1(strcmp(noninsbkps.chr1,chr) & noninsbkps.str1=='+' & noninsbkps.pos1<bkp1.pos1));
            if isempty(adjbkp_pos)
                adjbkp_pos=0; %set the nearest breakpoint to be the p-ter
            end
        else
            % for (+) breakpoints, find the (-) breakpoint to the right of
            % the last (+) breakpoint in the group
            adjbkp_pos=min(noninsbkps.pos1(strcmp(noninsbkps.chr1,chr) & noninsbkps.str1=='-' & noninsbkps.pos1>bkp2.pos1));
            if isempty(adjbkp_pos)
                adjbkp_pos=chrLengths(find(strcmp(ChrLabels,chr))); %set the nearest breakpoint to be q-ter
            end
        end
        
        PAR_loci=[PAR_loci;[bkp1(1,{'Sample','par_id','chr1','str1'}),dataset(no_bkp,bkp1.pos1,bkp2.pos1,adjbkp_pos,max(abs(adjbkp_pos-bkp2.pos1),abs(adjbkp_pos-bkp1.pos1)),...
            'VarNames',{'no_bkps','left_most_bkp','right_most_bkp','neighbor_bkp_pos','min_seg_length'})]];
    end
end

%% Adjacent parallel breakpoints with adjacent cis breakpoints
[Shared,IA,IB]=intersect(PAR(:,2:7),GAPS(:,1:6));
i=1;
A=[];B=[];
aid=1;bid=1;
while i<max(PAR.par_id)
    curr_par=PAR(PAR.par_id==i,:);
    if curr_par.str1(1)=='+' % adjacent (+/+) breakpoints 
        curr_bkp=curr_par(1,:);
        [~,~,cid]=intersect(curr_bkp(1,2:7),Shared);
        if ~isempty(cid)
            add=[GAPS(IB(cid)-1,:);curr_par(:,2:end)]; %look for adjacent (-) breakpoint to the left of the (+) breakpoint
            A=[A;dataset(repmat(aid,size(add(:,1))),'VarNames',{'cid'}),add];
            aid=aid+1;
        end
        i=i+1;
    else % adjacent (-/-) breakpoints
        curr_bkp=curr_par(end,:);
        i=i+1;
        next_par=PAR(PAR.par_id==i,:);
        next_bkp=next_par(1,:);
        if curr_bkp.Sample==next_bkp.Sample & strcmp(curr_bkp.chr1,next_bkp.chr1) & next_bkp.str1=='+' & next_bkp.pos1-curr_bkp.pos1<=dist_threshold 
            % adjacent (+/+) breakpoints near adjacent (-/-) breakpoints,
            % i.e., reciprocal parallel breakpoints
            B=[B;dataset(repmat(bid,size(curr_par(:,1))),'VarNames',{'cid'}),curr_par(:,2:end);dataset(repmat(bid,size(next_par(:,1))),'VarNames',{'cid'}),next_par(:,2:end)];
            bid=bid+1;
            i=i+1;
        else
            [~,~,cid]=intersect(curr_bkp(1,2:7),Shared);
            if ~isempty(cid)
                add=[curr_par(:,2:end);GAPS(IB(cid)+1,:)]; %look for adjacent (+) breakpoint to the right of the (-) breakpoint
                A=[A;dataset(repmat(aid,size(add(:,1))),'VarNames',{'cid'}),add];
                aid=aid+1;
            end
        end
    end
end

PAR_with_GAP=A; % adjacent parallel breakpoints with one adjacent cis breakpoint
PAR_with_PAR=B; % reciprocal adjacent parallel breakpoints

%save adj_with_cis_bkps.mat PAR*

%% Foldbacks and nested deletions
FB_with_GAP=[];
for i=1:max(PAR_with_GAP.cid)
    if length(unique(PAR_with_GAP.Idx(PAR_with_GAP.cid==i)))<length(PAR_with_GAP.Idx(PAR_with_GAP.cid==i))
        fbi(i)=1; % two or more adjacent breakpoints in a cluster form a junction (i.e., a foldback junction)
        FB_with_GAP=[FB_with_GAP;PAR_with_GAP(PAR_with_GAP.cid==i,:)];
    else
        fbi(i)=0;
    end
end

FB_with_FB=[];
FB_with_PAR=[];
NESTED_DEL=[];
DEL_plus_LR=[];
for i=1:max(PAR_with_PAR.cid)
    curr_cluster=PAR_with_PAR(PAR_with_PAR.cid==i,:);
    minus_cluster=curr_cluster(curr_cluster.str1=='-',:);
    plus_cluster=curr_cluster(curr_cluster.str1=='+',:);
    if length(unique(plus_cluster.Idx))<length(plus_cluster.Idx) | length(unique(minus_cluster.Idx))<length(minus_cluster.Idx)
        FB_with_PAR=[FB_with_PAR;curr_cluster]; % at least one foldback junction on either side
    end
    if length(unique(plus_cluster.Idx))<length(plus_cluster.Idx) & length(unique(minus_cluster.Idx))<length(minus_cluster.Idx)
        FB_with_FB=[FB_with_FB;curr_cluster]; % at least one foldback junction on each side
    end
    if length(intersect(minus_cluster.Idx,plus_cluster.Idx))>1
        NESTED_DEL=[NESTED_DEL;curr_cluster]; % at least two deletion junctions between the left and right AJB groups
    end
    if length(intersect(minus_cluster.Idx,plus_cluster.Idx))==1
        DEL_plus_LR=[DEL_plus_LR;curr_cluster]; % one deletion junction between the left and right AJB groups
    end
end

save FB_DEL.mat FB* NESTED_DEL DEL_plus_LR

%% Insertion analysis
Insertions=sortrows(dataset(INS.Sample(1:2:end),INS.chr1(1:2:end),INS.pos1(1:2:end),INS.pos1(2:2:end),'VarNames',{'Sample','chr','left','right'}),{'Sample','chr','left'});
%length(Insertions) 85,684

Left={};
Right={};
for si=1:length(Samples)
    junctions=Junctions(Junctions.Sample==si,:);
    junctions_2=junctions; 
        junctions_2.chr1=junctions.chr2;junctions_2.pos1=junctions.pos2;junctions_2.str1=junctions.str2;
        junctions_2.chr2=junctions.chr1;junctions_2.pos2=junctions.pos1;junctions_2.str2=junctions.str1;
        junctions_2.bkp1_type=junctions.bkp2_type;junctions_2.bkp2_type=junctions.bkp1_type;
    bkp=sortrows([junctions;junctions_2],{'chr1','pos1'}); % collect all breakpoints in a sample
    
    ins_bkps=Insertions(Insertions.Sample==si,:);
    chrs=unique(ins_bkps.chr,'stable');
    Left_adj_idx={};
    Right_adj_idx={};
    for ci=1:length(chrs)
        bkp_pos_chr=bkp(strcmp(bkp.chr1,chrs{ci}),:);
        ins_bkps_chr=ins_bkps(strcmp(ins_bkps.chr,chrs{ci}),:);
        left={};right={};
        for ki=1:length(ins_bkps_chr)
            left{ki}=bkp_pos_chr.Idx(bkp_pos_chr.pos1>=ins_bkps_chr.left(ki)-10000 & bkp_pos_chr.pos1<ins_bkps_chr.left(ki)); % breakpoint within 10kb to the left of the insertion
            right{ki}=bkp_pos_chr.Idx(bkp_pos_chr.pos1<=ins_bkps_chr.right(ki)+10000 & bkp_pos_chr.pos1>ins_bkps_chr.right(ki)); % breakpoint within 10kb to the right of the insertion
        end
        Left_adj_idx=[Left_adj_idx,left];
        Right_adj_idx=[Right_adj_idx,right];
    end
    Left=[Left,Left_adj_idx];
    Right=[Right,Right_adj_idx];
end
save Insertions.mat Insertions Left Right

adj=0;
for i=1:length(Insertions)
    if ~isempty(Left{i}) || ~isempty(Right{i})
        adj=adj+1;
    end
end
adj; % number of insertions with an adjacent breakpoint on either side

overlap=double(min(Insertions.right(1:end-1),Insertions.right(2:end))-Insertions.left(2:end)); % smaller insertion contained in a larger insertion ('nested')
overlap(diff(Insertions.Sample)~=0 | ~strcmp(Insertions.chr(1:end-1),Insertions.chr(2:end)))=nan;

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

save BKP_classification.mat A B C D E F G I bkp_density
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

