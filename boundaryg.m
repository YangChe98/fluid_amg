function [g1,g2]=boundaryg(x,y)
g1=function_u(x,y);
g2=function_v(x,y);

% if(x==1||x==0||y==0)
% %     g1=0;
% %     g2=0;
% g1=function_u(x,y);
% g2=function_v(x,y);
% elseif (y==1)
% %     g1=1;
% %     g2=0;
% g1=function_u(x,y);
% g2=function_v(x,y);
% else 
%     g1=function_u(x,y);
% g2=function_v(x,y);
% %     g1=0;
% %     g2=0;
% end
end