from pylab import*
from scipy.optimize import curve_fit

#define constant
D=2.
dx=1
dt=0.1
n=200

#define Gaussian function
def gaussian(X, mu, sig):
	return np.exp(-(X - mu)**2 / (2*sig**2))/(sqrt(2*pi*sig**2))

#define list to store time that the measurements take place
tnum=[500,1000,3000,5000,9000]
print "Please enter the time[s] when the measurements take place:"
tnum[0]= input('t1 = ')/dt
tnum[1]= input('t2 = ')/dt
tnum[2]= input('t3=  ')/dt
tnum[3]= input('t4=  ')/dt
tnum[4]= input('t5=  ')/dt
t=int(max(tnum))+10

#initlalize the density distribution
rho=zeros((t,2*n+1))
for i in range(n-5,n+6):
    rho[0][i]=1./10.

#define list to store 2n+1 x cordinate
x=linspace(-n*dx,n*dx,2*n+1)

#diffusion formula
for i in range(0,t-1):
	for j in range(1,2*n-1):
		rho[i+1,j]=rho[i,j]+D*dt*(rho[i,j+1]+rho[i,j-1]-2*rho[i,j])/(dx**2)

#define and initalize parametes for fit, then fit 
p=zeros((5,2))
coef=[[0,0],[0,0],[0,0],[0,0],[0,0]]
covg=[0,0,0,0,0]
fitfun=zeros((5,2*n+1))
for i in range(0,5):
	p[i]=[0,sqrt(2*D*tnum[i]*dt)]
	coef[i],covg[i] = curve_fit(gaussian,x,rho[tnum[i]],p[i])
	fitfun[i]=gaussian(x,coef[i][0],coef[i][1])
	print coef[i]

#make the plot for inital state
figure()
plot(x,rho[0])
grid()
xlabel('x axis')
ylabel('probability density')
savefig('initial.pdf')

#make the plot at later times
figure()
for i in range(0,5):
	plot(x,rho[tnum[i]],'r--')
	plot(x,fitfun[i],'c')
legend(['t1(num.)','t1(fit.)','t2(num.)','t2(fit.)','t3(num.)','t3(fit.)','t4(num.)','t4(fit.)','t5(num.)','t5(fit.)'], loc='best')
xlabel('x axis')
ylabel('probability density')
grid()
savefig('fit.pdf')

	
show()
