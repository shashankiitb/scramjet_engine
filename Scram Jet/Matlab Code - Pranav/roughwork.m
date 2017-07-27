
% clear global
% ans = 8/(1/2)*2;
% guess=[85*pi/180 asin(1/2)+pi/180]



guess=[85*pi/180 asin(1/8.81)+pi/180];    % ????Guess function is used to solve for numerical solution, needto define a interval.
beta=fzero(@ThetaBetaSolve,guess,[],6,8.81);
 