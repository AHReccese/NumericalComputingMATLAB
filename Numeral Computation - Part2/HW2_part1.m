First Problem
1.Newton Method
x0 = 0; %input('enter the approximate root: ');
syms f(x)
figure();
xlabel("|x(N)-x(N-1)|");
ylabel("f(x(N))");
f(x) = sin(x)*exp(x)+25*x+1;
df = diff(f,x);
index_N = 1;
while(true)
    x1 = vpa(x0 - (f(x0)/vpa(df(x0))));
    stem(abs(x1-x0),f(x1),'filled')
    hold on
    Trail_Members = sprintf('x(%d): %.10f we Have x(%d)-x(%d) that is:%.10f && f(x(%d)) is %.10f: ',index_N,x1,index_N,index_N-1,x1-x0,index_N,f(x1))
    if(f(x1) == 0 ||  abs(x1-x0)< 10^(-4)) % we use 10^(-4) as the HW says
        root = sprintf('root finded:%.10f \nNumber of goneStages: %d',x1,index_N)
        break;
    end
 index_N = index_N + 1; 
    x0 = x1;
end

2.Fixed point Method
x0 = 0; %input('enter the approximate root: ');
syms f(x)
figure();
xlabel("x(N)-x(N-1)");
ylabel("f(x(N))");
f(x) = (1-sin(x)*exp(x))/25;
index_N = 1;
while(true)
    x1 = vpa(f(x0));
    stem(abs(x1-x0),f(x1),'filled')
    hold on
    Trail_Members = sprintf('x(%d): %.10f we have x(%d)-x(%d) that is: %.10f && f(x(%d)) is %.10f',index_N,x1,index_N,index_N-1,x1-x0,index_N,f(x1))
    if(f(x1) == 0 || abs(x1-x0) < 10^(-4))  % we use 10^(-4) as the HW says
       root = sprintf('root finded: %.10f \nNumber of goneStages: %d',x1,index_N)
       break;
    end
index_N = index_N + 1;
x0 = x1;
end

3.secant Method
x0 = -1; %input('enter first guess of root');
x1 =  0; %input('enter second guess of root');
syms f(x)
figure();
xlabel("x(N)-x(N-1)");
ylabel("f(x(N))");
f(x) = sin(x)*exp(x)+25*x+1;
index_N = 1;
while(true)
    stem(abs(x1-x0),f(x1),'filled');
    hold on
    Trail_Members = sprintf('x(%d): %.10f we have x(%d)-x(%d) that is: %.10f && f(x(%d)) is %.10f',index_N,x1,index_N,index_N-1,x1-x0,index_N,f(x1))
    rootFinder = x1 - (f(x1)*(x1-x0))/(f(x1)-f(x0));
    if(abs(rootFinder - x1) < 10^(-4))
        Trail_Members = sprintf('x(%d): %.10f we have x(%d)-x(%d) that is: %.10f && f(x(%d)) is %.10f',index_N+1,rootFinder,index_N+1,index_N,rootFinder-x1,index_N+1,f(rootFinder))
        root = sprintf('root finded: %.10f \nNumber of goneStages: %d',rootFinder,index_N + 1)
        break
    else
        x0 = x1;
        x1 = rootFinder;
        index_N = index_N + 1;
    end
end

4.Bisection
a = -1; %input('enter a: ')
b =  0; %input('enter b: ')
error = 10^(-4); %input('enter error: ')
syms f(x)
figure();
xlabel("x(N)-x(N-1)");
ylabel("f(x(N))");
f(x) = sin(x)*exp(x)+25*x+1;
index_N = 1;
lastOne = 0;
while(true)
    c = (a+b)/2;
    wholeData = sprintf('a = %.10f && b = %.10f && c = %.10f',a,b,c)
    if(index_N >=2)
    Trail_Members = sprintf('x(%d): %.10f we have x(%d)-x(%d) that is: %.10f && f(x(%d)) is %.10f',index_N,c,index_N,index_N-1,c-lastOne,index_N,f(c))  
    stem(abs(c-lastOne),f(c),'filled')
    hold on
    end
    if((index_N >= 2) && (abs((c-lastOne))<error || f(c) == 0 ))
        root = sprintf('root finded: %.10f \nNumber of goneStages: %d',rootFinder,index_N)
        break;
    end
    if(f(c)*f(a) < 0)
        b = c;
    else
        a = c;  
    end 
    lastOne = c;
index_N = index_N + 1;    
end

5.Regula Falsi
a = -1; %input('enter a: ')
b =  0; %input('enter b: ')
error = 10^(-4); %input('enter error: ')
syms f(x)
figure();
xlabel("x(N)-x(N-1)");
ylabel("f(x(N))");
f(x) = sin(x)*exp(x)+25*x+1;
index_N = 1;
lastOne = 0;
nowOne = 0;
while(true)
    nowOne = (a*f(b)-b*f(a))/(f(b)-f(a));
    wholeData = sprintf('a = %.10f && b = %.10f && x(%d) = %.10f',a,b,index_N,nowOne)
    if(index_N >=2)
    Trail_Members = sprintf('x(%d): %.10f we have x(%d)-x(%d) that is: %.10f && f(x(%d)) is %.10f',index_N,nowOne,index_N,index_N-1,nowOne-lastOne,index_N,f(nowOne))  
    stem(abs(nowOne-lastOne),f(nowOne),'filled')
    hold on
    end
    if((index_N >= 2) && (abs((nowOne - lastOne)<error || f(nowOne) == 0 )))
        root = sprintf('root finded: %.10f \nNumber of goneStages: %d',rootFinder,index_N)
        break;
    end
    if(f(c)*f(a) < 0)
        b = c;
    else
        a = c;  
    end 
    lastOne = c;
index_N = index_N + 1;    
end

Second Problem
Solve The equation(finding the roots) via Newton Method
x0 = 0.5; %input('enter the approximate root: '); 0.5 was my quess
syms f(x)
figure();
xlabel("|x(N)-x(N-1)|");
ylabel("f(x(N))");
f(x) = x^3 - (3*x^2)*(2^(-x)) + 3*x*(4^(-x)) - 8^(-x);
df = diff(f,x);
index_N = 1;
while(true)
    x1 = vpa(x0 - (f(x0)/vpa(df(x0))));
    stem(abs(x1-x0),f(x1),'filled')
    hold on
    Trail_Members = sprintf('x(%d): %.10f we Have x(%d)-x(%d) that is:%.10f && f(x(%d)) is %.10f: ',index_N,x1,index_N,index_N-1,x1-x0,index_N,f(x1))
    if(f(x1) == 0 ||  abs(x1-x0)< 10^(-4)) % we use 10^(-4) as the HW says
        root = sprintf('root finded:%.10f \nNumber of goneStages: %d',x1,index_N)
        break;
    end
 index_N = index_N + 1; 
    x0 = x1;
end

Second problem(Question)
The answer of the asked Question:
why convergence is not quadratic 
lests take a look at the differential of the f(x) which is at the below :
 ( it is also a matlab output!)
as we know the newton Method converges to the root quadratically but as we read in 
the proof we thought that the differential of the main mathod(f(x)) is not zero but
in this example i will show you that by converging to the main root we will find out that 
the differential of the f(x) at x = root converges to zero
thus the newton method is not quadratically converged anymore (because we had thought that
diff(x) at x = root is not zero ( the root should be a simple root (not more!))
now i wanna show You that the x(n) seauences converges to the root that satisfies the equation 
diff(root (= lim x(n) when n --> infinity)  = 0
lets see ( i decreased the error from 10^-4 to 10^-8 in order to show that diff(x) converges to the zero when the 
x(n) sequences converges to the root (it is also an obvious fact in the graph)
x0 = 0.5; %input('enter the approximate root: '); 0.5 was my quess
syms f(x)
figure();
xlabel("|x(N)-x(N-1)|");
ylabel("f(x(N))");
f(x) = x^3 - (3*x^2)*(2^(-x)) + 3*x*(4^(-x)) - 8^(-x);
df = diff(f,x);
index_N = 1;
while(true)
    x1 = vpa(x0 - (f(x0)/vpa(df(x0))));
    stem(x1, diff(x1))
    hold on
    Trail_Members = sprintf('x(%d): %.6f we Have diff(x) that is:%.40f : ',index_N,x1,vpa(df(x1)))
    if(f(x1) == 0 ||  abs(x1-x0)< 10^(-8)) % we use 10^(-4) as the HW says
        root = sprintf('root finded:%.6f \nNumber of goneStages: %d',x1,index_N)
        break;
    end
 index_N = index_N + 1; 
    x0 = x1;
end





