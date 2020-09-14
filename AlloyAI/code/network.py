//============================DNN module=================================
Import random 
Import numpy as np 
Import json 

def rmse(x,y):
    z = []
    for i in range(len(x)):
        z.append((x[i] - y[i])**2)
    z_bar = np.mean(z)
    return np.sqrt(z_bar)
def mae(x,y):    
    z = []
    for i in range(len(x)):
        z.append(abs(x[i] - y[i]))
return np.mean(z)

class Network(object):
    def __init__(self,sizes):
        self.num_layers = len(sizes)
        self.sizes = sizes
        self.biases = [np.random.randn(b,1) for b in sizes[1:]]
        self.weights = [np.random.randn(x,y)/np.sqrt(y) for y, x in zip(sizes[:-1],sizes[1:])]
        #self.weights = [np.random.randn(x,y) for y, x in zip(sizes[:-1],sizes[1:])]
    def feedforward(self,a):
        icount = 0
        for b, w in zip(self.biases,self.weights):
            if icount < len(self.biases):
            #if icount < len(self.biases)-1:
                a = relu(np.dot(w,a)+b)
                icount += 1
            else:
                a = np.dot(w,a)+b
        return a
    def SGD(self,x_train,y_train,epochs,mini_batch_size,eta,lmbda=0.0,x_test=None,y_test=None):
        training_data = []
        test_data = []
        x_train = np.array(x_train)
        y_train = np.array(y_train)
        x_train = x_train.reshape(len(x_train),len(x_train[0]),1)
        y_train = y_train.reshape(len(y_train),1)
        for i in range(len(x_train)):
            training_data.append((x_train[i],y_train[i]))
        #for i in range(len(x_test)):
            #test_data.append((x_test[i],y_test[i])) 
        if test_data: n_test = len(test_data)
        n_train = len(training_data)
        for j in range(epochs):
            #eta_t = eta/pow(j+1,0.5)
            eta_t = eta
            random.shuffle(training_data)
            mini_batches = [training_data[k:k+mini_batch_size] for k in range(0,n_train,mini_batch_size)]
            for mini_batch in mini_batches:
                self.updata_mini_batch(mini_batch,eta_t,lmbda,len(training_data))
            if test_data:
                print("Epoch {0}: {1}".format(j,self.evaluate(test_data)))
            else:
                print("Epoch {0}: loss  {1} complete".format(j,self.evaluate(training_data)))
    def updata_mini_batch(self,mini_batch,eta,lmbda,n):
        momemtum = 0.9
        nabla_b = [np.zeros(b.shape) for b in self.biases]
        nabla_w = [np.zeros(w.shape) for w in self.weights]
        for x, y in mini_batch:
            delta_nabla_b, delta_nabla_w = self.backprop(x,y)
            nabla_b = [nb+dnb for nb,dnb in zip(nabla_b,delta_nabla_b)]
            nabla_w = [nw+dnw for nw,dnw in zip(nabla_w, delta_nabla_w)]
        self.weights = [(1-eta*(lmbda/n))*w-(eta/len(mini_batch))*nw for w,nw in zip(self.weights,nabla_w)]
        self.biases = [b-(eta/len(mini_batch))*nb for b,nb in zip(self.biases,nabla_b)]
    def backprop(self,x,y):
        nabla_b = [np.zeros(b.shape) for b in self.biases]
        nabla_w = [np.zeros(w.shape) for w in self.weights]
        activation = x
        activations = [x]
        zs = []
        for b,w in zip(self.biases,self.weights):
            z = np.dot(w,activation)+b
            zs.append(z)
            activation = relu(z)
            activations.append(activation)
        delta = self.cost_derivative(activations[-1],y)*relu_prime(zs[-1])
        nabla_b[-1] = delta
        nabla_w[-1] = np.dot(delta,activations[-2].transpose())
        for l in range(2,self.num_layers):
            z = zs[-1]
            sp = relu_prime(z)
            delta = np.dot(self.weights[-l+1].transpose(),delta)*sp
            nabla_b[-l] = delta
            nabla_w[-l] = np.dot(delta,activations[-l-1].transpose())
            return(nabla_b,nabla_w)
    def evaluate(self, data):
        y_pred = [self.feedforward(x) for (x,y) in data]
        y_orid = [y for (x,y) in data]
        result = rmse(y_pred,y_orid)
        return result
    def predict(self, xdata,ydata):
        xdata = np.array(xdata)
        ydata = np.array(ydata)
        xdata = xdata.reshape(len(xdata),len(xdata[0]),1)
        ydata = ydata.reshape(len(ydata),1)
        y_pred = np.array([self.feedforward(x) for x in xdata])
        y_orid = ydata
        y_pred = y_pred.reshape(len(y_pred))
        result = rmse(y_pred,y_orid)
        return result,y_pred,y_orid
    def cost_derivative(self,output_activations,y):
        return(output_activations-y)
    def save(self,filename):
        params = {'sizes':self.sizes,
                 'weights':[w.tolist() for w in self.weights],
                 'biases':[b.tolist() for b in self.biases]
                 }
        with open(filename,'w') as f:
            json.dump(params,f)
def load_model(filename):
    with open(filename,'r') as f:
        params = json.load(f)
    model = Network(params['sizes'])
    model.weights = [np.array(w) for w in params['weights']]
    model.biases = [np.array(b) for b in params['biases']]
    return model
def relu(z):
    #return (abs(z) + z) / 2.0
    z[z<0] =0
    return z
def relu_prime(z):
    z[z <= 0] = 0
    z[z > 0] = 1
return z
#==========================data collection=================================
def delete_repeat_row(a):
 new_a=[tuple(row) for row in a]
 s=set(new_a)
 return np.array(list(s))

def readbinaryunique(properties,s):
  #s=9, FeCr
  #s=10, FeAl
  #s=11, CrAl
    
  data= np.loadtxt('B:\\gcc\\bccCE'+str(s)+'.dat')
  data= np.delete(data,(4,5,6,7,8,9,10,11,13),1)
    
  new_data= delete_repeat_row(data)
  (l1,l2)=(np.shape(new_data))
    
  K=new_data[:,0].reshape((l1,1))
  G=new_data[:,1].reshape((l1,1))
  v=new_data[:,2].reshape((l1,1))
  E=new_data[:,3].reshape((l1,1))
  
  T=new_data[:,4].reshape((l1,1))
   if(s==10):
        xFe=new_data[:,5].reshape((l1,1))
        xCr=np.zeros((l1,1))
        xAl=np.ones((l1,1))-new_data[:,5].reshape((l1,1))
  elif(s==11):
        xFe=np.zeros((l1,1))
        xCr=new_data[:,5].reshape((l1,1))
        xAl=np.ones((l1,1))-new_data[:,5].reshape((l1,1))
  elif (s==9):
        xFe=new_data[:,5].reshape((l1,1))
        xCr=np.ones((l1,1))-new_data[:,5].reshape((l1,1))
        xAl=np.zeros((l1,1))

       
  x_train=np.hstack((xFe,xCr,xAl,T))
  if (properties in 'bulk modulus'):
      y_train=K
  elif (properties in 'shear modulus'):
      y_train=G
  elif (properties in 'Poisson ratio'):
        y_train=v
  elif (properties in 'Yong modulus'):
        y_train=E

  return x_train,y_train

def datatrainunique(properties):
  x1,y1=readbinaryunique(properties,9)
  x2,y2=readbinaryunique(properties,10)
  x3,y3=readbinaryunique(properties,11)
    
  x_train=np.vstack((x1,x2,x3))
  y_train=np.vstack((y1,y2,y3))
    
  return x_train,y_train

def readternaryunique(properties):
  data =np.loadtxt('B:\\gcc\\bccff3z.dat')
  data= np.delete(data,(5,6,7),1)
    
  data_test=delete_repeat_row(data)    
  (d1,d2)=(np.shape(data_test))

  composition=data_test[:,5:8].reshape((d1,3))

  T=data_test[:,4].reshape((d1,1))    
  K=data_test[:,0].reshape((d1,1))
  G=data_test[:,1].reshape((d1,1))
  v=data_test[:,2].reshape((d1,1))
  E=data_test[:,3].reshape((d1,1))
    
  x_test=np.hstack((composition,T))
    
  if (properties in 'bulk modulus'):        
        y_test=K
  elif (properties in 'shear modulus'):        
        y_test=G
  elif (properties in 'Poisson ratio'):
        y_test=v
  elif (properties in 'Yong modulus'):
        y_test=E
    
    for i in ['x','2','4','5','6']:
        data=np.loadtxt('B:\\gcc\\bccff3'+i+'.dat')
        data= np.delete(data,(5,6,7),1)
        data_test=delete_repeat_row(data)    
        (d1,d2)=(np.shape(data_test))
        
        composition=data_test[:,5:8].reshape((d1,3))

        T=data_test[:,4].reshape((d1,1))
        K=data_test[:,0].reshape((d1,1))
        G=data_test[:,1].reshape((d1,1))
        v=data_test[:,2].reshape((d1,1))
        E=data_test[:,3].reshape((d1,1))
        
        x2=np.hstack((composition,T))
        
        x_test=np.vstack((x_test,x2))
        if (properties in 'bulk modulus'):
            y_test=np.vstack((y_test,K))
        elif (properties in 'shear modulus'):
            y_test=np.vstack((y_test,G))
        elif (properties in 'Poisson ratio'):
            y_test=np.vstack((y_test,v))
        elif (properties in 'Yong modulus'):
            y_test=np.vstack((y_test,E))
    
    for i in ['a','b','c','d']:
        data=np.loadtxt('B:\\gcc\\bcczf3'+i+'.dat')
        data= np.delete(data,(17,18,19),1)
        data_test=delete_repeat_row(data)
        
        (d1,d2)=(np.shape(data_test))
        
        composition=data_test[:,17:20].reshape((d1,3))
        T=data_test[:,16].reshape((d1,1))
        K=data_test[:,12].reshape((d1,1))
        G=data_test[:,13].reshape((d1,1))
        v=data_test[:,14].reshape((d1,1))
        E=data_test[:,15].reshape((d1,1))
    
        x2=np.hstack((composition,T))
       
        x_test=np.vstack((x_test,x2))
        if (properties in 'bulk modulus'):
            y_test=np.vstack((y_test,K))
        elif (properties in 'shear modulus'):
            y_test=np.vstack((y_test,G))
        elif (properties in 'Poisson ratio'):
            y_test=np.vstack((y_test,v))
        elif (properties in 'Yong modulus'):
            y_test=np.vstack((y_test,E))    
        
    return x_test, y_test

def normailzation(x_train,x_test):
    scaler = StandardScaler()

    scaler.fit(x_train)
    x_train = scaler.transform(x_train)

    #scaler.fit(x_test)
    x_test = scaler.transform(x_test)
    
    return (x_train,x_test)

# train data and test data
def dataunique2(properties):
    x_train,y_train=datatrainunique(properties)
    x_test,y_test=readternaryunique(properties)    
    return x_train,y_train,x_test,y_test

def dataunique(properties):
    x_train,y_train=datatrainunique(properties)
    x_test,y_test=readternaryunique(properties)    
    #normalization
    scaler = StandardScaler()
    scaler.fit(x_train)
    x_train = scaler.transform(x_train)
    x_test = scaler.transform(x_test)
    return x_train,y_train,x_test,y_test


x_train,y_train,x_test3,y_test3=dataunique('bulk')
np.save("C:\\Users\\paopao\\Documents\\python2\\bulk\\x_binary_total.npy",x_train)
np.save("C:\\Users\\paopao\\Documents\\python2\\bulk\\y_binary_total.npy",y_train)
np.save("C:\\Users\\paopao\\Documents\\python2\\bulk\\x_ternary_total.npy",x_test3)
np.save("C:\\Users\\paopao\\Documents\\python2\\bulk\\y_ternary_total.npy",y_test3)

#=========================== training the model=============================
model = Network([4,11,9,10,1])
model.SGD(x_train2,y_train2,150,600,0.0002,0.0001)
a1,b1,c1=model.predict(x_test,y_test)
a2,b2,c2=model.predict(x_train,y_train)
x2=mae(b2,c2)
y2=rmse(b2,c2)
x1=mae(b1,c1)
y1=rmse(b1,c1)
print (x1,y1,x2,y2)

#=============================plot result ================================
from matplotlib import colors
from matplotlib import pyplot
from matplotlib.gridspec import GridSpec
import matplotlib.ticker as tk
import math
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties as FP

xmin = 10
xmax = 300
ymin = 10
ymax = 300
s = 1
a = 0.6
pyplot.rcParams['font.sans-serif'] = 'Times New Roman'
pyplot.rcParams['font.size'] = 18

fig = plt.figure(figsize=(8,4), dpi=72, facecolor="white") 
gs = GridSpec(1,2,width_ratios=[5,5])
gs.update(left=0.08,right=0.96,wspace=0.05,hspace=0.1)

efp = FP('Times New Roman', size=18)

#fig.text(0.41,0.00,r'$\rm{Experiment \ T_c}$',fontproperties=efp)
#fig.text(0.005,0.61,r'$\rm{Predicted \ T_c}$',rotation=90,fontproperties=efp)

fig.text(0.31,0.00,'Calculated $K_{VRH}$',fontproperties=efp)
fig.text(-0.03,0.61,'Predicted $K_{VRH}$',rotation=90,fontproperties=efp)
gs.update(left=0.06,right=1,wspace=0.1,hspace=0.1)

ax1 = plt.subplot(gs[0,0])
ax1.plot([xmin,xmax],[ymin,ymax],color='red',linewidth=1.0,linestyle='--')
ax1.scatter(c2.flatten(),b2.flatten(), edgecolor='b',c="blue", s=s, alpha=a)

plt.text(15, 250, 'Training',fontproperties=efp)
plt.xlim([xmin,xmax])
plt.ylim([ymin,ymax])
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

plt.grid(linestyle='--')
ax1.tick_params(labelsize=18)
labels = ax1.get_xticklabels() + ax1.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
ax1.xaxis.set_major_locator(tk.FixedLocator(np.arange(0,301,50)))

ax2 = plt.subplot(gs[0,1])

ax2.plot([xmin,xmax],[ymin,ymax],color='red',linewidth=1.0,linestyle='--')
ax2.scatter(c1.flatten(),b1.flatten(), edgecolor='b',c="blue", s=s, alpha=a)

plt.text(15, 250, 'Test',fontproperties=efp)
plt.xlim([xmin,xmax])
plt.ylim([ymin,ymax])
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.grid(linestyle='--')
ax2.yaxis.set_major_formatter(plt.NullFormatter())

ax2.tick_params(labelsize=18)
labels = ax2.get_xticklabels() + ax2.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]

[label.set_fontname('Times New Roman') for label in labels]
ax2.xaxis.set_major_locator(tk.FixedLocator(np.arange(0,301,50)))

#plt.style.use('classic')
plt.style.use('bmh')
plt.savefig("bulk_result.png",dpi=600,bbox_inches = 'tight')
plt.show()