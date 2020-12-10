I=zeros(201,1);
V=linspace(-1,1,201); 
for n=1:201
f=I(n)-10^(-11)*(exp((V(n)-I(n)*10)/0.026)-1.0);
df=1.0+10^(-10)/0.026*exp((V(n)-I(n)*10)/0.026);
I0(n)=abs(I(n)-f/df);

while abs((I0(n)-I(n))/(I0(n)+eps))>=0.001 
    I(n)=I0(n); 
    f=I(n)-10^(-11)*(exp((V(n)-I(n)*10)/0.026)-1.0); 
    df=1.0+10^(-10)/0.026*exp((V(n)-I(n)*10)/0.026); 
    I0(n)=abs(I(n)-f/df);

end
end
plot(V,I0);
axis([-1 1 -0.01 0.05]); 
xlabel('U/V'); 
ylabel('I/A');
