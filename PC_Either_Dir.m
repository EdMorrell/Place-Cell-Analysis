function [PC_Ver] = PC_Either_Dir(PC_Ver1,PC_Ver2)

%   ---    PC_Either_Dir
%              - Takes 2 place cell verify arrays (eg. left/right) and
%              returns verify array if cell a 'place-cell' in either
%              running direction

Pos1_PC = find(PC_Ver1);
Pos2_PC = find(PC_Ver2);
all_PC = [Pos1_PC;Pos2_PC];
all_PC = sortrows(all_PC);
pv = unique(all_PC);
PC_Ver = zeros(max(size(PC_Ver1,1),size(PC_Ver2,1)),1);
PC_Ver(pv) = 1;