clear
close all

B1(1) = 0.0000000000001;
df = 0;
dp = 0;
mode = 2;
mx = 0;
my = 0;
mz = 1;
gr(1) = 0;
gr(2) = 0;
t1 = 10;
    t2 = 1;

[outmx,outmy,outmz] = bloch([0 B1(1)],gr,[ 0 0.001],t1,t2,df,dp,mode,mx,my,mz)
signal(1) = (outmz(2));
for counter = 2 : 1500

    
    B1(counter) = B1(counter-1) * 1.02;
    [outmx,outmy,outmz] = bloch([ 0 0.001],[ 0 0],[0 B1(counter)],t1,t2,df,dp,mode,mx,my,mz)
    signal(counter) = (outmz(2));
    
end
plot((B1),signal)