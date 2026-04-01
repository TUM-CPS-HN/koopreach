function R_data = reach_LS(params,options)


% initialize cell array that stores the reachable sets
t = params.tStart:params.dt:params.tFinal;

steps = length(t)-1;
R_data = cell(steps+1,1);
R_data{1} = params.R0;
% loop over all reachablity steps
for i = 1:steps


    %reduce
    R_data{i} = reduce(R_data{i},'girard',100);
    % compute next reachable set

    % %%%%%%%%-------------------Data driven reachability-----------
    options.Uorig= params.U ;
    xStar = R_data{i}.center;
    uStar =options.Uorig.center;
    xStarMat = repmat(xStar,1,size(options.X_0T,2));
    uStarMat = repmat(uStar,1,size(options.U_full,2));
    oneMat = repmat([1],1,size(options.U_full,2));

    X1W_cen =  options.X_1T - options.Wmatzono.center;
    X1W = matZonotope(X1W_cen,options.Wmatzono.generator);
    % IAB = (options.X_1T )*pinv([oneMat; options.X_0T+(-1*xStarMat);options.U_full+-1*uStarMat]);
    IAB = (X1W )*pinv([oneMat; options.X_0T+(-1*xStarMat);options.U_full+-1*uStarMat]);

    V =  options.X_1T + -1*(IAB*[oneMat; options.X_0T+(-1*xStarMat);options.U_full+-1*uStarMat] );
    VInt = intervalMatrix(V);
    leftLimit = infimum(VInt);
    rightLimit = supremum(VInt);

    V_one= zonotope(interval(min(leftLimit')',max(rightLimit')'));

    R_data{i+1} = IAB*cartProd([1],cartProd(R_data{i}+(-1*xStar),options.Uorig+(-1*uStar))) +V_one+ options.W+options.Zeps ;

    

end

% create reachable set object
timePoint.time = num2cell(t(2:end)');
time = 2:3:length(t);

timePoint_data.set = R_data(time);
timePoint_data.time = num2cell(t(time)');

R_data = reachSet(timePoint_data);



end