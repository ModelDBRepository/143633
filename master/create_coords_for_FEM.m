
load Full_3d_VIM_locs_2009Feb.mat

nTC = 490;
nCER = 233;
nCTX = 233;
nRN = 233;
nTIN = 113;
NTOT = nTC+nCER+nCTX+nRN+nTIN;


% Arrange potentials from COMSOL by cell order
load Pop1XYZPhi_vim_0.2_UMFPACK.txt
PHI_CM1 = Pop1XYZPhi_vim_0(:,4);
PHI_NRN1 = cell(100,1);
Drop_ax1 = ones(5,100);
load Drop1_ind.mat

for nc=1:100
    
    PHI_NRN1{nc,1}(1:NTOT,1) = PHI_CM1((nc-1)*NTOT+1:(nc)*NTOT,1);
    
    miniDrop1 = find(Drop1((nc-1)*NTOT+1:(nc)*NTOT)==1);
    if isempty(miniDrop1)
    else
        if min(miniDrop1) <= nTC;
            Drop_ax1(1,nc) = 0;
        else
        end
        if ~isempty(find(miniDrop1 <= nTC + nCER & miniDrop1 > nTC));
            Drop_ax1(2,nc) = 0;
        else
        end
        if ~isempty(find(miniDrop1 <= nTC + nCER + nCTX & miniDrop1 > nTC + nCER));
            Drop_ax1(3,nc) = 0;
        else
        end
        if ~isempty(find(miniDrop1 <= nTC + nCER + nCTX + nRN & miniDrop1 > nTC + nCER + nCTX));
            Drop_ax1(4,nc) = 0;
        else
        end
        if max(miniDrop1) > nTC + nCER + nCTX + nRN;
            Drop_ax1(5,nc) = 0;
        else
        end
    end
end
save PHI_by_cell.mat PHI_NRN*

% Export Potentials for NEURON USE
load PHI_by_cell.mat

for n_cell = 1:100;
    f_name = (['phi_pop1_cell_' num2str(n_cell) '.dat']);
    fid = fopen(f_name,'wb');
    fwrite(fid,PHI_NRN1{n_cell,1},'double');
    fclose(fid);
end;

Drop_ax1 = Drop_ax1(1:500)
fid = fopen('Drop_ax1.dat','wb')
fwrite(fid,Drop_ax1,'double');
fclose(fid);






