permtest<-function(datamatrix, designmatrix, distance="euclid", nperms=1000, maxperms=1000000, designtype="random", savedist=FALSE, savedmat=FALSE)
{
 	designmatrix[,2]<-as.character(as.double(unlist(as.vector(t(designmatrix[,2]))))-1);
	designmatrix<-t(designmatrix);

#Remove rows from the datamatrix containing NaN's and NA's
	v<-vector();
	for(i in 1:nrow(datamatrix))
	{
		s<-sum(datamatrix[i,]);
		v[i]<-TRUE;
		if(is.nan(s) || is.null(s))
		{
			v[i]<-FALSE;
		}
	}
#t(as.matrix(datamatrix));
	datamatrix<-datamatrix[v,];
#check designmatrix for nan's/na's
	for(i in 1:nrow(designmatrix))
	{
		for(j in 1:ncol(designmatrix))
		{
			if(is.null(designmatrix[i,j]) || is.nan(unlist(designmatrix[i,j])))
			{
				#print("NaN's or NA's not allowed in designmatrix");
				return("NaN's or NA's not allowed in designmatrix");
			}
		}
	}
	
	if(nperms>maxperms)
	{
		nperms<-maxperms;
	}
	n2 <- sum(as.double(unlist(as.vector(t(designmatrix[2,])))));


	n1 <- ncol(designmatrix)-n2;

	n <- n1 + n2;
	dim <- nrow(datamatrix);
	ispaired<-FALSE;
	isblocked<-FALSE;
	if(!(designtype=="paired")&&!(designtype=="blocked")&&!(designtype=="random"))
	{
		#print("Unknown designtype");
		return("Unknown designtype");	
	}
	if(designtype=="paired")
	{
		ispaired<-TRUE;
	}
	if(designtype=="blocked")
	{
		isblocked<-TRUE;
	}
	colnames(designmatrix)<-designmatrix[1,];

	tdesignmatrix<-designmatrix;

	for(i in 1:ncol(datamatrix))
	{

		tdesignmatrix[1:nrow(tdesignmatrix),i]<-designmatrix[1:nrow(designmatrix),colnames(datamatrix)[i]];
	}
 	designmatrix<-tdesignmatrix; 
 	indexofsort<-sort(as.double(unlist(as.vector(t(designmatrix[2,])))),decreasing=FALSE,index.return=TRUE);
 	indexofsort<-indexofsort$ix;	
 	tdatamatrix<-datamatrix; 	 	
 	
 	for(i in 1:ncol(datamatrix))
	{
		tdatamatrix[1:nrow(tdatamatrix),i]<-datamatrix[1:nrow(datamatrix),indexofsort[i]]

		tdesignmatrix[1:nrow(tdesignmatrix),i]<-designmatrix[1:nrow(designmatrix),indexofsort[i]]
	}
 	designmatrix<-tdesignmatrix; 
	datamatrix<-tdatamatrix;
	colnames(designmatrix)<-designmatrix[1,];
	colnames(datamatrix)<-designmatrix[1,];
#check for missing data items
	if((sum(as.integer(colnames(datamatrix)==colnames(designmatrix))) != ncol(datamatrix)))
	{
		#print("Designmatrix is incomplete or contains duplicates")
		return("Designmatrix is incomplete or contains duplicates")
	}

#sort groups 0 and 1 by pairing	
 	if(ispaired)
 	{
	 	tdesignmatrix<-designmatrix; 
		tdatamatrix<-datamatrix;
		for(i in 1:n1)
		{
			tdesignmatrix[1:nrow(tdesignmatrix),i]<-designmatrix[1:nrow(designmatrix),i];
			tdatamatrix[1:nrow(tdatamatrix),i]<-datamatrix[1:nrow(datamatrix),i];	
			tdesignmatrix[1:nrow(tdesignmatrix),(i+n1)]<-designmatrix[1:nrow(designmatrix),designmatrix[[3,i]]]
			tdatamatrix[1:nrow(tdatamatrix),(i+n1)]<-datamatrix[1:nrow(datamatrix),designmatrix[[3,i]]];		
		}
		colnames(tdesignmatrix)<-tdesignmatrix[1,];
		colnames(tdatamatrix)<-tdesignmatrix[1,];
		designmatrix<-tdesignmatrix;
		datamatrix<-tdatamatrix;
 	}

#sort by block and by group within block 	
	if(isblocked)
	{
	 	indexofsort<-sort(as.character(designmatrix[3,]),index.return=TRUE);
	 	indexofsort<-indexofsort$ix;
	 	tdatamatrix<-datamatrix; 
	 	tdesignmatrix<-designmatrix
	 	for(i in 1:ncol(datamatrix))
		{
			tdatamatrix[1:nrow(tdatamatrix),i]<-datamatrix[1:nrow(datamatrix),indexofsort[i]]
			tdesignmatrix[1:nrow(tdesignmatrix),i]<-designmatrix[1:nrow(designmatrix),indexofsort[i]]
		}	 		 	
		colnames(tdesignmatrix)<-tdesignmatrix[1,];
		colnames(tdatamatrix)<-tdesignmatrix[1,];
	 	designmatrix<-tdesignmatrix;
	 	datamatrix<-tdatamatrix;
		blknames<-unique(as.character(designmatrix[3,]));
		nblks<-sum(nchar(nchar(nchar(nchar(nchar(nchar(blknames)))))));
		blkctl<-matrix(data=0,nrow=5+n,ncol=nblks);
		currindx<-0;
		oldindx<-1;
		for(i in 1:nblks)
		{
			currindx<-currindx+1;
			currblkname<-designmatrix[3,currindx];
			for(j in currindx:ncol(designmatrix))
			{
				if((as.character(designmatrix[3,j]) != currblkname)||(j==ncol(designmatrix)))
				{
					oldindx<-currindx;
					currindx<-(j-1);
					if(j==ncol(designmatrix))
					{
						currindx<-j;
					}
					locn1<-sum(as.double(designmatrix[2,oldindx:currindx]));
					locn2<-currindx-(oldindx-1)-locn1; 
					blkctl[1,i]<-oldindx;
					blkctl[2,i]<-currindx;
					blkctl[3,i]<-locn1;
					blkctl[4,i]<-locn2;
					blkctl[5,i]<-locn1*locn1-locn1;
					blkctl[6,i]<-locn2*locn2-locn2;
					blkctl[7,i]<-locn1*locn2;	
					for(k in 1:locn1)
					{
						blkctl[(7+k),i]<-1;
					}
					break
				}
			}
			
		}
		blkctl[5,]<-sum(blkctl[5,]);
		blkctl[6,]<-sum(blkctl[6,]);
		blkctl[7,]<-sum(blkctl[7,]);
	}
		

#compute the initial distance matrix
	if(distance=="logr")
	{
		dmat<-matrix(nrow=ncol(datamatrix),ncol=ncol(datamatrix));
		dmat<- -(log(cor(datamatrix)));
		for(i in 1:ncol(datamatrix))
		{
			for(j in 1:ncol(datamatrix))
			{
				if(is.nan(dmat[i,j]))
				{
#print("negative logs detected, must use distance:euclid")
					distance<-"euclid";
					return("negative logs detected, must use distance:euclid");
				}
			}
		}
	}
	if(distance=="euclid")
	{	
		n<-ncol(datamatrix);
		dmat<-matrix(nrow=n,ncol=n);

		for (i in 1:n)
		{
			for (j in 1:n)
			{
				dmat[i,j]<-((t(unlist(datamatrix[1:nrow(datamatrix),i])-unlist(datamatrix[1:nrow(datamatrix),j]))))%*%((unlist(datamatrix[1:nrow(datamatrix),i])-unlist(datamatrix[1:nrow(datamatrix),j])))/nrow(datamatrix);
		    	}		
		}
	}
		
	if(ispaired)
	{
		for(i in 1:n2)
		{
			dmat[(n1+i),i]<-0;
			dmat[i,(n1+i)]<-0;
		}
	}





	
	if(isblocked)
	{
		tdmat<-matrix(data=0,nrow=n,ncol=n);
		for(i in 1:nblks)
		{
			for(j in blkctl[1,i]:blkctl[2,i])
			{
				for(k in blkctl[1,i]:blkctl[2,i])
				{
					tdmat[j,k]<-dmat[j,k];
				}
			}
		}
		dmat<-tdmat;
	}
	
	
	
########################################	
#      	cell means and test statistics
########################################
	if(isblocked)
	{
		runn11<-0;
		runn22<-0;
		runn12<-0;
		for(i in 1:nblks)
		{
			n1<-blkctl[3,i];
			n2<-blkctl[4,i];
			n<-n1+n2
			i1<-matrix(nrow=(n1+n2),ncol=2);
			i1[1:n1,1]<-blkctl[8:(7+n1),i];
			if(n2 != 0)
			{
				i1[(n1+1):(n1+n2),1]<-0;
			}
			i1[,2]<-!i1[,1];
			tdmat<-dmat[blkctl[1,i]:blkctl[2,i],blkctl[1,i]:blkctl[2,i]];
			sums<-(t(i1)%*%tdmat%*%i1);
			runn11<-runn11+sums[1,1];
			runn22<-runn22+sums[2,2];
			runn12<-runn12+sums[1,2];			
		}
#divide by global
			if(runn11!=0)
			{
				mn11<-runn11/blkctl[5,1];
			}else
			{
				mn11<-runn11;
			}
			if(runn22!=0)
			{
				mn22<-runn22/blkctl[6,1];
			}else
			{
				mn22<-runn22;
			}
			if(runn12!=0)
			{
				mn12<-runn12/blkctl[7,1];
			}else
			{
				mn12<-runn12;
			}
		rdiff<-mn11-mn22;
		adiff<-abs(mn22-mn11);
		sum<-mn12-(mn11+mn22)/2;		
		omnibus<-mn12-mn11;	

	}else
	{

		n<-n1+n2
		n11<-n1*n1-n1;
		n22<-n2*n2-n2;
		n12<-n1*n2;
		
		if(n11==0)
		{
			n11=1;
		}
		if(n22==0)
		{
			n22=1;
		}
		if(n12==0)
		{
			n12=1;
		}

		if(ispaired)
		{	
			n12<-((n1*n2)-n2);
		}
		i2 <- rnorm((n1+n2)*(n1+n2),0,0);
		dim(i2) <- c((n1+n2),(n1+n2));
		i3 <- rnorm((n1+n2)*(n1+n2),1,0);
		dim(i3) <- c((n1+n2),(n1+n2));
		i1 <- i3[1:(n1+n2),1:2];            


		i1[(n1+1):(n1+n2),1] <- i2[(n1+1):(n1+n2),1];           
		i1[1:n1,2] <- i2[1:n1,2];  

		sums<-(t(i1)%*%dmat%*%i1);

		mn11<-sums[1,1]/n11;
		mn22<-sums[2,2]/n22;
		mn12<-sums[1,2]/n12;

		rdiff<-mn11-mn22;
		adiff<-abs(mn22-mn11);
		sum<-mn12-(mn11+mn22)/2;		
		omnibus<-mn12-mn11;
	}
	
	countrdiff<-0;
	countadiff<-0;
	countsum<-0;
	countomnibus<-0;
	distrdiff<-list();
	distadiff<-list();
	distsum<-list();
	distomnibus<-list();
	
############################################perms###################################

	if(isblocked)
	{
		for (z in 1:nperms)
		{
			prunn11<-0;
			prunn22<-0;
			prunn12<-0;
			for(i in 1:nblks)
			{
				n1<-blkctl[3,i];
				n2<-blkctl[4,i];
				n<-n1+n2
				n11<-n1*n1-n1;
				n22<-n2*n2-n2;
				n12<-n1*n2;

				i1<-matrix(nrow=(n1+n2),ncol=2);
				i1[1:n1,1]<-blkctl[8:(7+n1),i];
				i1[(n1+1):(n1+n2),1]<-0;
				i1[1:(n1+n2),1]<-sample(i1[1:(n1+n2),1]);
				i1[,2]<-!i1[,1];

				tdmat<-dmat[blkctl[1,i]:blkctl[2,i],blkctl[1,i]:blkctl[2,i]];
				psums<-(t(i1)%*%tdmat%*%i1);
				prunn11<-prunn11+psums[1,1];
				prunn22<-prunn22+psums[2,2];
				prunn12<-prunn12+psums[1,2];	
			}
#divide by global
			if(prunn11!=0)
			{
				pmn11<-prunn11/blkctl[5,1];
			}else
			{
				pmn11<-prunn11;
			}
			if(prunn22!=0)
			{
				pmn22<-prunn22/blkctl[6,1];
			}else
			{
				pmn22<-prunn22;
			}
			if(prunn12!=0)
			{
				pmn12<-prunn12/blkctl[7,1];
			}else
			{
				pmn12<-prunn12;
			}

			prdiff<-pmn11-pmn22;
			padiff<-abs(pmn22-pmn11);
			psum<-pmn12-(pmn11+pmn22)/2;		
			pomnibus<-pmn12-pmn11;
			
			countsum<-countsum+(psum>=sum);
			countadiff<-countadiff+(padiff>=adiff);
			if(rdiff <= 0)
			{
				countrdiff<-countrdiff+(prdiff<=rdiff);
			}else
			{
				countrdiff<-countrdiff+(prdiff>=rdiff);			
			}
			countomnibus<-countomnibus+(pomnibus>=omnibus);    
	
			if(savedist)
			{
				distrdiff[[z]]<-prdiff;
				distadiff[[z]]<-padiff;
				distsum[[z]]<-psum;
				distomnibus[[z]]<-pomnibus;
			}
		}
	}
	if(ispaired)
	{
		if(n1<13)
		{
			nperms<-(2^n1);
			pil<-list();
			for(i in 1:((2^n1)-1))
			{
				pit<-matrix(nrow=1,ncol=n1);
				k<-i;
				for(j in (n1):1)
				{
					kt<- k-(2^(j-1));
					if(kt>=0)
					{
						pit[j]<-TRUE;
						k<-kt;
					}
					if(kt<0)
					{	
						pit[j]<-FALSE;
					}
				}
				pil[[i]]<-t(t(as.integer(pit)));
			}
			pil[[nperms]]<-t(t(as.integer(!(pil[[((2^n1)-1)]]))));
			for (i in 1:nperms)
			{
				pi<-matrix(nrow=(n1+n2),ncol=2);
				pit<-pil[[i]];
				pi[1:n1,1]<-pit[1:n1,1];
				pi[(n1+1):(n1+n2),1]<-!pit[1:n1,1];
				pi[,2]<-!pi[,1];
				psums<-t(pi)%*%dmat%*%pi;

			    	pmn11<-psums[1,1]/n11;
			    	pmn22<-psums[2,2]/n22;
			    	pmn12<-psums[1,2]/n12;

				prdiff<-pmn11-pmn22;
				padiff<-abs(pmn22-pmn11);

				psum<-pmn12-(pmn11+pmn22)/2;		
				pomnibus<-pmn12-pmn11;

				countsum<-countsum+(psum>=sum);
				countadiff<-countadiff+(padiff>=adiff);
				if(rdiff <= 0)
				{
					countrdiff<-countrdiff+(prdiff<=rdiff);
				}else
				{
					countrdiff<-countrdiff+(prdiff>=rdiff);			
				}
				countomnibus<-countomnibus+(pomnibus>=omnibus);    
		
				if(savedist)
				{
					distrdiff[[i]]<-prdiff;
					distadiff[[i]]<-padiff;
					distsum[[i]]<-psum;
					distomnibus[[i]]<-pomnibus;
				}
			}
		}
		if(n1>=13)
		{
#print("Total cols exceed 13, switching from enumerations to random sampling")
			nperms<-maxperms;
			for (i in 1:maxperms)
			{
				pit<-matrix(nrow=1,ncol=n1);
				for(j in 1:n1)
				{
					pit[[j]]<-sample(matrix(data=list(1,0),nrow=2,ncol=1))[[1]];
				}
				pit<-t(pit);
				pi<-matrix(nrow=(n1+n2),ncol=2);
				pi[1:n1,1]<-pit[1:n1,1];
				pi[(n1+1):(n1+n2),1]<-!pit[1:n1,1];
				pi[,2]<-!pi[,1];
				psums<-t(pi)%*%dmat%*%pi;

			    	pmn11<-psums[1,1]/n11;
			    	pmn22<-psums[2,2]/n22;
			    	pmn12<-psums[1,2]/n12;

				prdiff<-pmn11-pmn22;
				padiff<-abs(pmn22-pmn11);
				psum<-pmn12-(pmn11+pmn22)/2;		
				pomnibus<-pmn12-pmn11;
	    
				countsum<-countsum+(psum>=sum);
				countadiff<-countadiff+(padiff>=adiff);
				if(rdiff <= 0)
				{
					countrdiff<-countrdiff+(prdiff<=rdiff);
				}else
				{
					countrdiff<-countrdiff+(prdiff>=rdiff);			
				}
				countomnibus<-countomnibus+(pomnibus>=omnibus);    
		
				if(savedist)
				{	
					distrdiff[[i]]<-prdiff;
					distadiff[[i]]<-padiff;
					distsum[[i]]<-psum;
					distomnibus[[i]]<-pomnibus;
				}
			}
		}
	}
	if((!ispaired)&&(!isblocked))
	{	
		for (i in 1:nperms)
		{
			pi<-matrix(nrow=(n1+n2),ncol=2);
			pi[,1]<-sample(i1[,1]);;
			pi[,2]<-as.integer(!pi[,1]);	    
			psums<-t(pi)%*%dmat%*%pi;

		    	pmn11<-psums[1,1]/n11;
		    	pmn22<-psums[2,2]/n22;
		    	pmn12<-psums[1,2]/n12;
		    	
			prdiff<-pmn11-pmn22;
			padiff<-abs(pmn22-pmn11);
			psum<-pmn12-(pmn11+pmn22)/2;		
			pomnibus<-pmn12-pmn11;
			countsum<-countsum+(psum>=sum);
			countadiff<-countadiff+(padiff>=adiff);
			if(rdiff <= 0)
			{
				countrdiff<-countrdiff+(prdiff<=rdiff);
			}else
			{
				countrdiff<-countrdiff+(prdiff>=rdiff);			
			}
			countomnibus<-countomnibus+(pomnibus>=omnibus);    
			
			if(savedist)
			{
				distrdiff[[i]]<-prdiff;
				distadiff[[i]]<-padiff;
				distsum[[i]]<-psum;
				distomnibus[[i]]<-pomnibus;
			}
		}
	}
	sumpvalue<-countsum/(nperms);
	adiffpvalue<-countadiff/(nperms);
	rdiffpvalue<-countrdiff/(nperms);
	omnibuspvalue<-countomnibus/(nperms);
	out<-list();
	DESCRIPTION<-c(
	"variability test(one sided)  :  ",
	"variability test(two sided)  :  ",
	"location test                :  ",
	"equivalence test             :  ",
	"mean within group 1 distance :  ",
	"mean within group 2 distance :  ",
	"mean between group distance  :  ",
	"number of permutations       :  ",
	"distance type                :  ",
	"design type                  :  ");
	COMPUTED<-c("d11-d22","abs(d11-d22)","d12-((d22+d11)/2)","d12-d11","d11","d22","d12",nperms,distance,designtype);
	STAT<-c(rdiff,adiff,sum,omnibus,mn11,mn22,mn12,NA,NA,NA)
	COUNT<-c(countrdiff,countadiff,countsum,countomnibus,NA,NA,NA,NA,NA,NA)
	PVALUE<-c(rdiffpvalue,adiffpvalue,sumpvalue,omnibuspvalue,NA,NA,NA,NA,NA,NA)
	out<-data.frame(DESCRIPTION,COMPUTED,STAT,COUNT,PVALUE);
	if(savedist || savedmat)
	{
		if(savedist)
		{
			out2<-list(v1=distrdiff,v2=distadiff,loc=distsum,equiv=distomnibus);
			return(list(table=out,distributions=out2));
		}
		if(savedmat)
		{
			return(list(table=out,distancematrix=dmat));
		}
		if(savedist && savedmat)
		{
			out2<-list(v1=distrdiff,v2=distadiff,loc=distsum,equiv=distomnibus);
			return(list(table=out,distributions=out2,distancematrix=dmat));
		}
	}
	return(out);
}