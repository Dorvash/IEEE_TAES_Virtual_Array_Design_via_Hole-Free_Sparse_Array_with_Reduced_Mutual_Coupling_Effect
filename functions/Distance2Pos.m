%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
%
% coding by Dr. Jianyan LIU 
%    with the School of Information and Electronics, Beijing Institute of Technology, China,
%    also with the School of Electrical and Electronic Engineering, Nanyang Technological University, Singapore 
%    (e-mail: jianyanliu@qq.com )
% Fun: conver array distance list to array position list 
% 
function [ C  ] = Distance2Pos( d1 )
N = length(d1) + 1 ;
C = zeros(1, N);
for i = 2 : N 
    C(i) = C(i-1) + d1(i-1);
end