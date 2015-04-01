from pylab import*
from scipy.optimize import curve_fit

#define constant
D=2.
dx=1
dt=0.1
n=200
t=10000

#define Gaussian function
def gaussian(X, mu, sig):
	return np.exp(-(X - mu)**2 / (2*sig**2))/(sqrt(2*pi*sig**2))

#define list to store probalbility density for 2n+1 points
rho=zeros((t,2*n+1))

#define list to store 2n+1 x cordinate
x=linspace(-n*dx,n*dx,2*n+1)

#define list to store time that the measurements take place
tnum=[500,1000,3000,5000,9000]

#initlalize the density distribution
for i in range(n-5,n+6):
	rho[0][i]=1./10.

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
legend(['500s(num.)','500s(fit.)','1000s(num.)','1000s(fit.)','3000s(num.)','3000s(fit.)','5000s(num.)','5000s(fit.)','9000s(num.)','9000s(fit.)'], loc='best')
xlabel('x axis')
ylabel('probability density')
grid()
savefig('fit.pdf')

	
show()
