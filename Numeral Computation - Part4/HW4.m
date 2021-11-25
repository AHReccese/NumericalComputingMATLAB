First Question
matrixDimension = 2;
for index = 1:10
    process = sprintf("stage%d : epsilon = 10 ^(-%d)",index,2*index)
    epsilon = ((10)^((-2)*index))
    AwithB = [[epsilon,1;1,1],[1 + epsilon;2]]

% Guass Method
AwithB = guass(AwithB,matrixDimension);
% now we have a diagonal matrix

CalResult = zeros(2,1);
    for i = 1:2
    CalResult(i,1) = (AwithB(i,matrixDimension + 1)/AwithB(i,i));   
    end
CalResult   
MainResult = [1;1]
stem(CalResult(1,1),CalResult(2,1))
hold on
end
Question answer : 
as we see in the plot when epsilon decreases, first we feel that the main answer(result) doesn’t have a noticeable change but when
we reach into epsilon = 10 ^ (-16) then suddenly we see that X1 jumps into 2.22 and then in the next step we see that X1 jumps into 
111.02  again in the next step we see that X1 jumps into 11102 so we can conclude that when epsilon decreases too much then 
 X1 increases suddenly but sharply, however X2  doesn’t change too much.
and the main reason of this fact is that when we don’t pivot the matrix then we will have much cancellation error (while
epsilon decreases) and as epsilon decreases  cancellation error during division and multiplication increases much as the fact 
we saw in the plot.
if we pivot the matrix then we will see that the result matrix doesnt change!
to see that just change the matrix AwithB = [[epsilon,1;1,1],[1 + epsilon;2] to
AwithB = [[1,1;epsilon,1],[2;1+epsilon]] then you will see that even at the last step the resutl is: 
transpose[1,1].

Second Question
function's implementation is at the bottom of the code
% simple example
A = [2,-1,-2;-4,6,3;-4,-2,8];
determinant(A,3);   % function's implementation is at the bottom of the code

Third Question
matrixDimension = 2;
figure
% samples which are NOT around of epsilonMahchine
% some examples of epsilon from : 0.01 until 0.00000000000000000001 (10 ^ (-16)) 
for index = 1:6
    process = sprintf("stage%d : epsilon = 10 ^(-%d)",index,index);
    epsilon = vpa((10)^((-1)*index));
    AwithB = [[1,1 + epsilon;1 - epsilon,1],[1 + (1 + epsilon)*epsilon;1]];
    justA = [1,1 + epsilon;1 - epsilon,1];
% Gauss Method
AwithB = guass(AwithB,matrixDimension);      
% now we have a diagonal matrix
CalResult = vpa(zeros(2,1));
ResultComponentsReleativeError = vpa(zeros(2,1));
mainResult = [1;epsilon];
    for i = 1:2
    CalResult(i,1) = vpa(AwithB(i,matrixDimension + 1)/AwithB(i,i));  
    ResultComponentsReleativeError(i,1) =vpa(sqrt(abs((CalResult(i,1))^2 - (mainResult(i,1))^2 )));
    end
mainResult
CalResult
ResultComponentsReleativeError
conditionNumberOfMat = conNum(justA)
stem(CalResult(1,1),CalResult(2,1),'filled');
hold on
stem(ResultComponentsReleativeError(1,1),ResultComponentsReleativeError(2,1));
hold on
end
% now some epsilons around sqrt(Machine Epsilon)!!!
% i find out that My epsilon machine is around 5.42e-20;
% i try the calculation for sqrt(MachEps) , sqrt(MachEps)/100 , sqrt(MachEps)/10 , 10*sqrt(MachEps) , 100*sqrt(MachEps) 
% lets go
figure
epsilonMachine = vpa(5.42e-20);
for index = -2:2
    process = sprintf("stage %d : epsilon = sqrt(MachineEpsilon) * 10 ^(%d)",index,index)
    epsilon = sqrt(epsilonMachine)*(10 ^ index);
    AwithB = [[1,1 + epsilon;1 - epsilon,1],[1 + (1 + epsilon)*epsilon;1]]; 
    justA = AwithB(:,1:2);
% Gauss Method
AwithB = guass(AwithB,matrixDimension);      
% now we have a diagonal matrix
CalResult = vpa(zeros(2,1));
ResultComponentsReleativeError = vpa(zeros(2,1));
mainResult = [1;epsilon];
    for i = 1:2
    CalResult(i,1) = vpa(AwithB(i,matrixDimension + 1)/AwithB(i,i));  
    ResultComponentsReleativeError(i,1) =vpa(sqrt(abs((CalResult(i,1))^2 - (mainResult(i,1))^2 )));
    end
mainResult
CalResult
ResultComponentsReleativeError
conditionNumberOfMat = conNum(justA)
stem(CalResult(1,1),CalResult(2,1),'filled');
hold on
stem(ResultComponentsReleativeError(1,1),ResultComponentsReleativeError(2,1));
hold on
end
Asked Questions Answers:
1. How accurately is each component determined?What conclusions can you draw from this experiments? 
Answer: when the epsilon is too  more than sqrt(machine epsilon) as we detemined in the first part we got that in this time
condition Number of the Matrix is quite little and relative error in each component of the solution is less to
but as we get nearer to sqrt(machine sqrt) as we did at the last examples of first part and mainly at the second part
we got that relative error in each component of the solution becomes more and more and mainly the condition Number
of matrix increases sharply too... So as we get nearer to sqrt(machine epsilon) Not only condition Number of  matrix 
increases sharply but the also relative error in each component of the solution increases too.
and we can also see this fact on th graphs!

Fourth Question
%%%%%%%%%% Fucntions's implementations are at the end of Code ( Functions's part! )
% examples : 
% first part a:
A1 = [0.641,0.242;0.321,0.121];
A2 = [10,-7,0;-3,2,6;5,-1,5];
[CalculatedNorm2_INVA1,CalculatedNormInf_INVA1] = findConditionNumber(A1,2)
MATLAB_NORM2_OF_invA1 = vpa(norm(inv(A1),2));
MATLAB_NORM_INF_OF_invA1= vpa(norm(inv(A1),'inf'));
[CalculatedNorm2_INVA2,CalculatedNormInf_INVA2] = findConditionNumber(A2,3);
MATLAB_NORM2_OF_invA2 = vpa(norm(inv(A2),2));
MATLAB_NORM_INF_OF_invA2 = vpa(norm(inv(A2),'inf'));
% second part b
% choose random V vectors!
yVecsA2 = zeros(3,5); 
yVecsA2(:,1) = [2;4;0];
yVecsA2(:,2) = [0.5;0.8;0];
yVecsA2(:,3) = [1.6;0.3;0];
yVecsA2(:,4) = [0.34;1.4;0];
yVecsA2(:,5) = [0.2;0.4;0];
yVecsA1 = zeros(2,5); 
yVecsA1(:,1) = [2;4];
yVecsA1(:,2) = [0.5;0.8];
yVecsA1(:,3) = [1.6;0.3];
yVecsA1(:,4) = [0.34;1.4];
yVecsA1(:,5) = [0.2;0.4];
ZA1vectors = zeros(2,5)
ZA2vectors = zeros(3,5)
for j = 1:5
    ZA1vectors(:,j) = mldivide(A1,yVecsA1(:,j))
    ZA2vectors(:,j) = mldivide(A2,yVecsA2(:,j))
end
CondNumsA1_2 = zeros(1,5);
CondNumsA1Inf = zeros(1,5);
CondNumsA2_2 = zeros(1,5);
CondNumsA2Inf = zeros(1,5);
for i = 1:5
    CondNumsA1_2(1,i) = vpa(norm(ZA1vectors(:,i),2)/norm(yVecsA1(:,i),2));
    CondNumsA1Inf(1,i) = vpa(norm(ZA1vectors(:,i),'inf')/norm(yVecsA1(:,i),'inf'));
    CondNumsA2_2(1,i) = vpa(norm(ZA2vectors(:,i),2)/norm(yVecsA2(:,i),2));
    CondNumsA2Inf(1,i) = vpa(norm(ZA2vectors(:,i),'inf')/norm(yVecsA2(:,i),'inf'));
end
MATLAB_NORM2_OF_invA2
MATLAB_NORM_INF_OF_invA2
MATLAB_NORM2_OF_invA1
MATLAB_NORM_INF_OF_invA1
Cond2NumA2 = vpa(max(CondNumsA1_2))
CondinfNumA2 = vpa(max(CondNumsA2Inf))
CondinfNumA1 = vpa(max(CondNumsA1Inf))
Cond2NumA1 = vpa(max(CondNumsA2_2))
as we see the first method's precision is more than second method's precision which we
used random matrixs in and the first method's calculation is more accurate than second one!

Functions's part
function det = determinant(mat,dimension)  %%% LU factorizaiotn:second Question determinant calculater function
lower = zeros(dimension,dimension);
upper = zeros(dimension,dimension);
for i = 1:dimension
    lower(i,i) = 1;
end

for i = 1:dimension 
        % Upper Triangular 
        for k = i:dimension 
  
            % Summation of L(i,j) * U(j,k) 
         sum = 0; 
            for j = 1:i-1
                sum = sum + (lower(i,j) * upper(j,k));
            end
            % Evaluating U(i,k) 
            upper(i,k) = mat(i,k) - sum;
        end
        % Lower Triangular 
        for k = i:dimension
            if (i == k) 
                lower(i,i) = 1; % Diagonal as 1     
                else  
                % Summation of L(k, j) * U(j, i) 
                sum = 0; 
                    for j = 1:i-1 
                    sum = sum + (lower(k,j) * upper(j,i)); 
                    end
                % Evaluating L(k, i) 
                lower(k,i) = (mat(k,i) - sum) / upper(i,i); 
            end
        end
end

upper  % display lower matrix
lower  % diplay  upper matrix
detupper = 1; % lower matrix's determinant 
detlower = 1; % upper matrix's determinant
for i = 1:dimension
    detlower = detlower*lower(i,i);
    detupper = detupper*upper(i,i);
end
det  = detlower * detupper  % as the fact that if A = B * C we have: det(A) = det(B)*det(C).
end

%%%%%%%%%%%%  ConditionNumber

function conditionNumber = conNum(mat)
infNormofMat = norm(mat,'inf');
infNormofMatInverse = norm(inv(mat),'inf');
conditionNumber = infNormofMat * infNormofMatInverse;
end

%%%%%%%%%%%%  Guass Method

function guassMat = guass(AwithB,matrixDimension)

for k=1:matrixDimension-1  % this fat FOR makes all the elements of the A matrix which are under the main diameter ZERO
        for i=k+1:matrixDimension
            m=AwithB(i,k)/AwithB(k,k);
            for j=k:matrixDimension+1
            AwithB(i,j)=AwithB(i,j)-m*AwithB(k,j);
            end
        end
end

   for k=1:matrixDimension-1   % this fat FOR makes all the elements of the A matrix which are upper the main diameter ZERO
        for i=1:k
            m=vpa(AwithB(i,k+1)/AwithB(k+1,k+1));
            for j=k+1:matrixDimension+1
            AwithB(i,j)=AwithB(i,j)-m*AwithB(i+1,j);
            end
        end
   end
guassMat = AwithB;   
end

%%%%%%%%%%    solving the equation: Transpose(U) * Y = C; 

function AnswerVector = findTriangleAnswer(UTranspose,C,dimension)

AnswerVector = zeros(dimension,1);
sum =  0;
for i = 1:dimension 
    for j = 1:(i-1)
        sum = sum + UTranspose(i,j) * AnswerVector(j,1);
    end
AnswerVector(i,1) = vpa(C(i,1)-sum)/UTranspose(i,i); 
sum = 0;
end

end





%%%%%%%%    findConditionNumber    
function [NORM2,NORMINF] = findConditionNumber(A1,A1dimension)

[~,U1] = lu(A1)
coloumnIndexNum = 2 ^ A1dimension;
% create all C vectors
Cmat = vpa(zeros(A1dimension,coloumnIndexNum)); % the last line is for vector's condition number
Vmat = vpa(zeros(A1dimension + 1,coloumnIndexNum));
for i = 1:A1dimension
    for j = 1:coloumnIndexNum
        Cmat(i,j) = 1;
    end
end
for i = 1:A1dimension
    for j = 1:(coloumnIndexNum)
        Cmat(i,j) = Cmat(i,j) * ((-1)^(floor(j/(coloumnIndexNum/(2^(i))))));
    end
end
Cmat

for i = 1:coloumnIndexNum
Vmat(1:A1dimension,i) = findTriangleAnswer(transpose(U1),Cmat(:,i),A1dimension); 
end
% now lets find the vector which is more in magnitude
for i = 1:coloumnIndexNum
TempSum = 0;     
    for j = 1:A1dimension
    TempSum = TempSum + Vmat(j,i)^2; 
    end
    Vmat(A1dimension + 1,i) = vpa(sqrt(TempSum));  
end
Vmat
maxMagnitude = max(Vmat(A1dimension + 1,:));
for i = 1:coloumnIndexNum
    if( Vmat(A1dimension + 1,i) == maxMagnitude )
      % mainV = Vmat(1:2,i);
        mainC = Cmat(1:A1dimension,i);
        break
    end
end
Ymat = vpa(mldivide(transpose(A1),mainC));
foundZ = mldivide(A1,Ymat);
NORM2 = norm(foundZ,2)/norm(Ymat,2);
NORMINF = norm(foundZ,'inf')/norm(Ymat,'inf');
end

