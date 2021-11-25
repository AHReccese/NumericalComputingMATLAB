Third Question
I used Newton Method as a zero finder to determine the elapsed time required to fall a distance of 1 Km. 
g = 9.8065;
k = 0.00341;
x0 = 50; %input('enter the approximate root: ');
syms f(x)
figure();
xlabel("|x(N)-x(N-1)|");
ylabel("f(x(N))");
f(x) = log((cosh(x*sqrt(g*k))))/k - 1000;
df = diff(f,x);
index_N = 1;
while(true)
    x1 = vpa(x0 - (f(x0)/vpa(df(x0))));
    stem(abs(x1-x0),f(x1),'filled')
    hold on
    Trail_Members = sprintf('x(%d): %.19f we Have x(%d)-x(%d) that is:%.19f && f(x(%d)) is %.19f: ',index_N,x1,index_N,index_N-1,x1-x0,index_N,f(x1))
    if(f(x1) == 0 ||  abs(x1-x0)< 10^(-4)) % we use 10^(-4) as the HW says
        root = sprintf('root finded:%0.19f \nNumber of goneStages: %d',x1,index_N)
        break;
    end
 index_N = index_N + 1; 
    x0 = x1;
end

Forth Question

Main solution of the Matrix via Matlab own Functions!
A = vpa(0);  % in order to make the variables of matrix in vpa mode
B = vpa(0);  % in order to make the variables of matrix in vpa mode
for i=1:5
   for j=1:5
        A(i,j)=vpa(abs((1/3)^(abs(i-j))));
    end
end
for i=1:5
    B(i,1)=vpa((1/7)^i);
end
A % in order to show A
B % in order to show B
% The result
answer = mldivide(A,B)
Gauss-Jordan
AwithB = [A,B]   % in order to use Gauss-Jordan iteration
GaussAnswer=vpa(0); % matrix of soloution
n=5; % dimension of the A matrix
for i=1:n-1
    for j=i+1:n
        if abs(AwithB(j,i))>abs(AwithB(i,i))
            T=AwithB(j,:);
            AwithB(j,:)=AwithB(i,:);
            AwithB(i,:)=T;
        end
    end
end
disp('After pivoting');
disp(AwithB);
for k=1:n-1
    for i=k+1:n
        m=AwithB(i,k)/AwithB(k,k);
        for j=k:n+1
            AwithB(i,j)=AwithB(i,j)-m*AwithB(k,j);
        end
    end
end
disp('Triangularize Form  ');
disp(AwithB);
            
if AwithB(n,n)==0
    disp('No unique solution');
end
    GaussAnswer(n)=AwithB(n,n+1)/AwithB(n,n);
    for j=n-1:-1:1
        sum=0;
        for i=1:n-j
            sum=sum+AwithB(j,n+1-i)*GaussAnswer(n+1-i);
        end
        GaussAnswer(j)=(AwithB(j,n+1)-sum)/AwithB(j,j);
    end 
    GaussAnswer
   % the error
   for i = 1:5
       sprintf('error of the Gauss jordan(%d element) from the main result of matlab is %.32f:\n and the gauss %d element is: %.32f and the main matlab %d result is %.32f ',i,vpa(abs(answer(i,1)-GaussAnswer(1,i))),i,GaussAnswer(1,i),i,answer(i,1))
   end

jacobi
Ajacobi = vpa(0);  % in order to make the variables of matrix in vpa mode
Bjacobi = vpa(0);  % in order to make the variables of matrix in vpa mode
for i=1:5
   for j=1:5
        Ajacobi(i,j)=vpa(abs((1/3)^(abs(i-j))));
    end
end
for i=1:5
    Bjacobi(i,1)=vpa((1/7)^i);
end
N=100; % 100 iteration
diagonal = diag(diag(Ajacobi)); % strip out the diagonal
diag_deleted = Ajacobi - diagonal; % delete the diagonal	
sol = vpa(0); % initial guess of zero
temp = sol;
    for i = 1:N
    % computing the matrix inverse    
    temp = diagonal \ (Bjacobi - diag_deleted * sol);
    end
sol = temp;
% the error
disp('')
    for i = 1:5
       sprintf('error of the Jacobi(%d element) from the main result of matlab is %.32f:\n and the jacobi %d element is: %.32f and the main matlab %d result is %.32f ',i,vpa(abs(answer(i,1)-sol(i,1))),i,sol(i,1),i,answer(i,1))
    end

Gauss-Seidel
ASeidel = vpa(0);  % in order to make the variables of matrix in vpa mode
BSeidel = vpa(0);  % in order to make the variables of matrix in vpa mode
for i=1:5
   for j=1:5
        ASeidel(i,j)=vpa(abs((1/3)^(abs(i-j))));
    end
end
for i=1:5
    BSeidel(i,1)=vpa((1/7)^i);
end
xSeidel=[vpa(0);vpa(0);vpa(0);vpa(0);vpa(0)];
SeidelAnswer=GaussSeidel(ASeidel,bb,xSeidel,100)
for i=1:5
   sprintf('error of Gauss-Seidel(%d element) from the main result of matlab is %.32f:\n and the seidel %d element is: %.32f and the main matlab %d result is %.32f',i,vpa(abs(answer(i,1)-SeidelAnswer(i,1))),i,SeidelAnswer(i,1),i,answer(i,1)) 
end


function y=GaussSeidel(AA,bb,xx,NumIters)
% Runs the Gauss-Seidel method for solving Ax=b, starting with x and
% running a maximum of NumIters iterations.
%
% The matrix A should be diagonally dominant, and in particular,  it should
% not have any diagonal elements that are zero (a division by zero error
% will be produced).
%
% The output y will be the whole sequence of outputs instead of the final
% value (if x is in R^n, then y will be n x NumIters
D=diag(AA);
AA=AA-diag(D);
D=vpa(1./D);  %We need the inverses
n=length(xx);
xx=xx(:); %Make sure x is a column vector
SeidelIter=zeros(n,NumIters);
for j=1:NumIters
for k=1:n
xx(k)=vpa((bb(k)-AA(k,:)*xx)*D(k));
end
SeidelIter(:,j)=xx
end
y=xx;
end


