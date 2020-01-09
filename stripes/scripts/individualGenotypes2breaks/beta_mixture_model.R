#
# 
# Author: patel
###############################################################################
args <- commandArgs(trailingOnly = TRUE)
sliding_window_outcome<-args[1]
outfile<-args[2]
cutoffs<-function(p1,p2,alpha,beta)
{
  return(beta(alpha+p1,beta+p2)/beta(alpha,beta));
}

beta_prob<-function(y,alpha_v,beta_v)
{
  alpha_v[alpha_v==0]=100000000
  beta_v[beta_v==0] = 100000000
  
  ##print(alpha_v)
  ##print(beta_v)
  return ((y^(alpha_v-1)*(1-y)^(beta_v-1))/beta(alpha_v, beta_v))
}

beta_prob_2<-function(y,alpha_v,beta_v,p)
{
  ## if(alpha_v <0 | beta_v<0)
  ## {
  ##     return (-10000000000000);
  ## }
  loglike<-(log(y^(alpha_v-1)*(1-y)^(beta_v-1)/beta(alpha_v, beta_v)))+log(p)
  return (loglike)
}
mix_beta_model_log_function<-function(x,data,z,p_)
{
  alpha<-x[1:3];
  beta<-x[4:6];
  #total_sum<-c();
  local_sum<-sapply(data,beta_prob_2,alpha=x[1:3],beta=x[4:6],p=p_[])	
  ll<-sum(t(local_sum)*z)*-1
  return(ll)
}


mix_beta_model<-function(x,data,z,p_)
{
  ##cat("now running it with ",x," ",p,"\n",file="data1.txt", sep="\t", append=TRUE)
  alpha<-x[1:3];
  beta<-x[4:6];
  ##print(x)
  midterm<-c();
  for(j in 1:3)
  {	   
    local_sum<-0;
    local_sum<-log(beta_prob(data,alpha[j],beta[j]))
    ##print (log(beta_prob(data,alpha[j],beta[j])))
    local_sum<-(local_sum+log(p_[j]));
    midterm<-append(midterm,local_sum);
  }
  ll<-sum(z*midterm)*(-1);
  #ll<-sumi*(-1);
  return (ll);	
}


dominator_nominator_func<-function(y,p,alpha,beta)
{
  return((p*beta_prob(y,alpha,beta)))
}


e_step<-function(x,data,p)
{
  alpha<-x[1:3];
  beta<-x[4:6];
  new_z<-c()
  dominator<-sapply(data,dominator_nominator_func,p=p[],alpha=x[1:3],beta=x[4:6])
  new_z<-apply(dominator, 1, function(x) x / colSums(dominator))
  return (new_z)
}

prob<-function(z)
{
  return(colSums(z)/dim(z)[1])
}

m_step<-function(x,data,z,p)
{
  v<-tryCatch({
    suppressWarnings(nlm(f=mix_beta_model,p=x,data=data,p_=p,z=z,print.level=0))
    },
    error = function(cond){
      return(NA)
  },
  warning= function(cond){
    ##print(cond)
    return(NULL)
  },
  finally={}			
  )
  return (v);
}
#####

main_f<-function(filename)
{
  
  d<-read.table(filename);
  
  #sample down to 50%
  k<-as.data.frame(table(d[,3]))
  k[,2]<-k[,2]*0.5
  data<-rep(as.numeric(levels(k[,1]))[k[,1]],k[,2])
  data <- data + 100
  data <- data/200
  data[data==0] = 0.00001
  data[data==1] = 1-0.00001
  #d[,3]<-d[,3]+100
  #d[,3]<-d[,3]/200
  #d[d[,3]==0,3]<-0.00001 #beta function is not defined at 0 and 1 
  #d[d[,3]==1,3]<-1-0.00001
  
  #data<-d[,3];
  z<-c();
  
  
  z<-matrix(0,length(data),3)
  for(i in 1:length(data))
  {
    if(data[i]>0.75)
    {
      z[i,]<-c(0.001,0.001,0.998);
    }else if(data[i]<0.25)
    {
      z[i,]<-c(0.998,0.001,0.001);
    }else{
      z[i,]<-c(0.001,0.998,0.001);
    }
  }
  
  big_counter<-0;
  x<-c(1,1,1,1,1,1);
  while(TRUE & big_counter < 20){

    counter<-0;
    v<-c();
    diff_log=1000000;
    v_pre = NA
    while(TRUE & counter < 100){
      ##cat("counter: ",counter," \n");
      #z<-model.matrix(~0 + as.factor(z));
      p<-prob(z);
      ##cat("p ",p,"\n");
      #m-step
      #v<-nlm(f=mix_beta_model,p=x,data=data,p_=p,z=z)
      #v<-nlm(f=mix_beta_model,p=x,data=data,p_=p,z=z,#print.level=2)
      v<-m_step(x,data,z,p)
      ##print(x)
      ##print(z)
      ##print(p)
      ##cat("oof\n")
      ##print(v)
      ##cat("v-pre")
      ##print(v_pre)
      
      if(anyNA(v)) ## if v is NA
      {		         
        if(anyNA(v_pre)){ ## ... and this is the first iteration
          #cat("v is na, and v_pre is still empty, breaking \n");
          #print(v)
          break ## break, since there will be no good v_pre to fall back to
        }
        v<-v_pre; ## else, set v to older v (v_pre)
        #cat("v is NA, falling back to v_pre, break \n ")
        break;
      }
      else
      {
        v_pre<-v; # if v 
        #cat("v was good, setting v_pre to v, continue \n ")
	#print(v)
      }
      if(length(v$estimate[v$estimate>50]) & counter>10) 
      {
        #cat("there's at least one estimate above 50, break.");
	#print(v$estimate);
	break;
      }
      
      #v<-nlm(f=mix_beta_model_log_function,p=x,data=data,p_=p,z=z)
      #e-step
      #z<-e_step(v$par,data,p)
      z<-e_step(v$estimate,data,p) 
      x<-v$estimate;
      #cat("computing new z and x \n")
      if(abs(v$minimum-diff_log)<0.2){
        #cat("v doesnt change significantly between iterations, break\n");
	break;
      	}
      diff_log<-v$minimum;
      counter<-counter+1;
      ##cat(v)
      #produce the histogram 
      #cat("end loop nr: \n")
      #print(counter)
      }

      if(!anyNA(v)){
        #cat("v is not NA, break \n");break
	}
      else{
	   x <- runif(6,0.99,1.01); # randomise the starting probabilities, see if it converges
      	   big_counter = big_counter+1
      	   #cat("didnt converge, using starting_probabilities\n");
      	   #print(x)
      	  }
    } #big_loop
    
    #return (v)
    #while
    #output
    if(anyNA(v)){ # if the model cant converge, do this: 
      upper_bound <- NULL
      lower_bound <- NULL
    }
    else {
      x<-v$estimate
      #rownames(z)<-d[,3];
      #t<-z[z[,1]<z[,2] & z[,2]>z[,3],]
      s<-seq(0.001,0.99,by=0.001)
      k<-matrix(0,4,length(s))
      k[1,]<-s
      k[2,]<-dbeta(s,x[1],x[1+3])
      k[3,]<-dbeta(s,x[2],x[2+3])
      k[4,]<-dbeta(s,x[3],x[3+3])
      k<-t(k);
      #LRN: To avoid confusion I changed the name of the matrix t to t1
      #LRN: To avoid any floating point errors I added a small threshold to be reached before the beta distributions could be considered different
      #LRN: Added the check that the vector t1 is not empty
      float.err1 = 1e-12; t1 <- k[ k[,3] - k[,2] > float.err1 & k[,3] - k[,4] > float.err1,] #LRN
      ## Added by tpayen, case we found in our dataset
      ## If there is no recombination in the dataset (scaffold in this case)
      ## t1 is unidimentional, we check for it
      if(!is.null(dim(t1))){
      #cat("t1 is okay\n")
        upper_bound<-max(t1[,1])
        upper_bound<-upper_bound*200-100
        lower_bound<-min(t1[,1])
        lower_bound<-lower_bound*200-100
      }
      else{
      #cat("t1 is null\n")
        upper_bound <- NULL
        lower_bound <- NULL
      }	
    } # else for if the model converges
    
    #now we have to take care, if values are to extreme.......
    #based on experience
    if(length(lower_bound)==0 || lower_bound >0)
    {
      #cat("modifying lower bound due to extreme values...\n")
      lower_bound<--25;
      
    }else if(lower_bound < -75)
    {
      #cat("modifying lower bound due to extreme values...\n")
      lower_bound<- -50;
    }
    
    if(length(upper_bound)==0 || upper_bound> 90)
    {
      #cat("modifying upper bound due to extreme values...\n")
      upper_bound<-50;
    }else if(upper_bound < 10)
    {
      #cat("modifying upper bound due to extreme values...\n")
      upper_bound<- 25;
    }
    
    border_matrix<-matrix(0,2,1)
    border_matrix[1,]<-lower_bound
    border_matrix[2,]<-upper_bound
    write.table( border_matrix, sep="\n", file=outfile, row.names=F,col.names=F)
    
    
  }
  
  
  main_f(sliding_window_outcome)
  
  
  
  
  
  ## s<-seq(0.001,0.99,by=0.001)
  ## xx<-v$estimate;
  ## plot(s,dbeta(s,xx[3],xx[3+3]),col="blue")
  ## points(s,dbeta(s,xx[1],xx[1+3]),col="red")
  ## points(s,dbeta(s,xx[2],xx[2+3]),col="green")
  ## points(s,dbeta(s,xx[3],xx[3+3]),col="blue")