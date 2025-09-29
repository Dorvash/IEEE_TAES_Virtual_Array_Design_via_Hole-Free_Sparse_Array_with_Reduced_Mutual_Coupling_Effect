function [ pos,  eleNum,  virApr ] = AugNestedII1Gen( N , L1, L2)
%   ref: Jianyan LIU,...etc  
%     "Augmented Nested Arrays with Enhanced DOF and Reduced Mutual Coupling" 
%     submitted to IEEE Transactions on Signal Processing , 2017.
%
%coding by Dr. Jianyan LIU 
%    with the School of Information and Electronics, Beijing Institute of Technology, China,
%    also with the School of Electrical and Electronic Engineering, Nanyang Technological University, Singapore 
%    (e-mail:jianyanliu@qq.com)
%
%fun: Augmented nested array II
%
%input arg: 
%     N      @ elememt number 
%     L1     @ the element number of right sub-ULA
%     L2     @ the element number of left sub-ULA
%output arg: 
%     pos    @ element position of physical NLA 
%     eleNum @ elememt number 
%     virApr @ virtual aperture of co-array
%
r1 =  round((N-4)/6);
r2 = r1 ;
r3 = r2 ; 
N1 = 4 * r1 + 3 ;
N2 = N -  N1;
if  (nargin == 3)
   r1 = L1;
   r2 = L2;
   r3 =  max(L1,L2);
   N1 =  r1 +  r2 + 2 * r3 + 3 ;
   N2 = N - N1 ;
end
d1=[ones(1,r1)  (r1+r2+2) ones(1,r3)*(r1+r2+2) ones(1,N2)*(2*r1+2*r2+3) ones(1,r3)*(r1+r2+1) (r1+1)   ones(1,r2)  ];
pos = Distance2Pos(d1);
pos = unique(sort(pos ));
eleNum = length(pos);
virApr =  max(pos) - min(pos) ;
