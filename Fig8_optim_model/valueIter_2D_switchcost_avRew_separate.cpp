// Compilar con -mcmodel=medium para que compile cuando la memoria es grande
// g++ valueIter_2D_switchcost_avRew_separate.cpp -mcmodel=medium , en cluster
// icpc runs but gives weird results (??)

#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define MAXDV           (5) //3
#define DVDELTA 		(0.1)
#define NSTEPS 		    (40) // 60


// #define SIGMA 			(0.1414) //corresponde a var=1 en 1 sec;depende de dt
// float deltatime = 0.025; //en sec. 0.01
float deltatime = 0.05; //en sec. 0.01
//#define DELTATIME       (0.020);//en secs
//#define KAPPA           (13) //9
#define KAPPA           (5.5) //9

#define NDV 			(2* int(MAXDV/DVDELTA) + 1) // number of mu values
#define NBEL  		    (2*NDV*NSTEPS*NDV*NSTEPS) //2 due to the tracking of last action
const int Nactions = 10; // 8 terminal actions; 2 listens


bool force_alternation = 1; //0,1
#define UNSIGNED (1) //0,1
 


#if UNSIGNED==0
    #define NCOH 12
#else
    #define NCOH 6
#endif

float coherences[NCOH] = {};

#if (UNSIGNED==0)
    float metaCoherences[NCOH] = {-0.64, -0.512, -0.384, -0.256, -0.128, -0.0000001,0.0000001, 0.128, 0.256, 0.384, 0.512, 0.64};
#else
    float metaCoherences[NCOH] = {0.0000001, 0.128, 0.256, 0.384, 0.512, 0.64};
#endif



//#define NCOH (12)
//float coherences[NCOH] = {};
//float metaCoherences[NCOH] = {-0.64, -0.512, -0.384, -0.256, -0.128, -0.0000001,0.0000001, 0.128, 0.256, 0.384, 0.512, 0.64};


float p_prior_coh[NCOH] = {};


float SIGMA = sqrt(deltatime);

//#define SIGMA 			(0.1414) //corresponde a var=1 en 1 sec;depende de dt

int Beliefs[NBEL][5]; //mu1,step1, mu2, step2, iLast
float V[NBEL];
float Vprev[NBEL];
float Vdelta[NBEL];
float T[NBEL];
int estadoNan[NBEL];

float Q[NBEL][Nactions];
int BestAction[NBEL];
int updateOrder[NBEL];

int   TlistenInd[NBEL][NDV][2]; //el dos es por las dos posibles acciones de listen
float TlistenVal[NDV][NSTEPS][NDV];

int   TStartInd;
float TStartVal;

//float Tstart[NBEL];
float Rewards[NBEL][Nactions];

float dv[NDV];

//float sigmaDV[NSTEPS];

float rewardListen      = 0.0;//0.2
//float rewardSuccess     = 1.0;//20
//float rewardFailure     = 0.0;
float rewardTimesUp     = -1.0;
// float Rdiff = 0.0;
// float Reasy = 1.0;
// float Rhard = 1.0;

float Rdiff = 1.0;
float Reasy = 0.0;
float Rhard = 0.0;

float DummyPenalty = -5.0;



        
        
//float Tintertrial   = 0.5;// en segundos
float Tintertrial   = 1.0;// en segundos
float Tpenalty      = 1.0;// en segundos
float TnonDecision  = 0.33;//en segundos
// float TswitchDim    = 0.1; //en segundos
float TswitchDim    = 0.0; //en segundos
//float TswitchDim    = 0.2; //en segundos
float rho;//reward por paso

float gama 		= 1.0;
float maxError 	= 0.00001;//0.000001;

float qact[Nactions];


const double pi = 3.14159265;

double normalpdf(double x,double mu, double sigma)
{  //You will change this function if you want another type of distribution
    
    return exp( -1 * (x - mu) * (x - mu) / (2 * sigma * sigma)) / (sigma * sqrt(2 * pi));
}


// normalize vector to sum to 1
void normalizeToOne(float *V,int n)
{
    float suma = 0.0;
    for (int i=0;i<n;i++){
        suma += V[i];
    }
    
    
    for (int i=0;i<n;i++){
        V[i] *= 1.0/suma;
        //std::cout << V[i] << '\n';
    }
    
//    if (V[0]!=V[0])
//        std::cout << suma << '\t';

}



void calc_pcoh(float *p_coh, float DV, int timestep, float coherences[])
{
    float media,desvio;
    
    int i;

    for (i=0;i<NCOH;i++){
        
        media = coherences[i]*((float) timestep);
        desvio = (float) SIGMA*((sqrt( (float) timestep)));
        
        if (timestep>0)
            p_coh[i] = (float) p_prior_coh[i] * normalpdf(DV,media,desvio);
        else
            p_coh[i] = (float) p_prior_coh[i];
    }
    
    normalizeToOne(p_coh,NCOH);
    
//    if (p_coh[0]!=p_coh[0]) // nan -> return prior
//    {
//        for (i=0;i<NCOH;i++){
//            p_coh[i] = (float) p_prior_coh[i];
//        }
//        normalizeToOne(p_coh,NCOH);
//    }
    
}


float calc_pplus(float DV, int timestep, float coherences[])
{
    float p_positive;
    float p_coh[NCOH];
    int i;
    
    calc_pcoh(p_coh,DV,timestep,coherences);
    
    // see positives
    p_positive = 0;
    for (i=0;i<NCOH;i++){
        if (coherences[i]>0)
            p_positive += p_coh[i];
        else if (coherences[i]==0)
            p_positive += p_coh[i]/2.0;
    }
    return p_positive;
}




void saveToFile() {
	std::ofstream myfile;
	myfile.open ("output.txt");
	myfile << "mu1\t step1\t mu2\t step2\t last_query\t value\t bestAction\n" ;
    
	for (int i = 0; i < NBEL; i++) {
		myfile << dv[Beliefs[i][0]] << '\t';
		myfile << Beliefs[i][1] << '\t';
        myfile << dv[Beliefs[i][2]] << '\t';
        myfile << Beliefs[i][3] << '\t';
        myfile << Beliefs[i][4] << '\t';
        
		myfile << V[i] << '\t';
		myfile << BestAction[i] << '\n';
    }
    
    myfile.close();
    
}

void savePars() {
    // creates a matlab ready file
    std::ofstream f;
    f.open("params.m");
    f << "pars.MAXDV = " << MAXDV <<";\n";
    f << "pars.DVDELTA = " << DVDELTA <<";\n";
    f << "pars.NSTEPS = " << NSTEPS <<";\n";
    f << "pars.SIGMA = " << SIGMA <<";\n";
    f << "pars.KAPPA = " << KAPPA <<";\n";
    f << "pars.deltatime = " << deltatime <<";\n";
    f << "pars.nonDecisionTime = " << TnonDecision <<";\n";
    f << "pars.switchDimTime = " << TswitchDim <<";\n";
    f << "pars.Tintertrial = " << Tintertrial <<";\n";
    f << "pars.Tpenalty = " << Tpenalty <<";\n";

    f << "pars.unsigned_bool = " << UNSIGNED <<";\n";
    f << "pars.force_alternation_bool = " << force_alternation <<";\n";
    
    
//    f << "pars.task_flag = " << task_flag <<";\n";

    f << "pars.EV_STEP = [" << coherences[0];
    for (int i=1;i<NCOH;i++){
        f << ",";
        f << coherences[i];
    }
    f << "];\n";
    
    f << "pars.COH = [" << coherences[0]/(deltatime*KAPPA);
    for (int i=1;i<NCOH;i++){
        f << ",";
        f << coherences[i]/(deltatime*KAPPA);
    }
    f << "];\n";
    
    f.close();
}

void saveQ(){
    std::ofstream myfile;
    myfile.open ("Qba.txt");
    for (int i = 0; i < NBEL; i++) {
        for (int j = 0;j<Nactions;j++) {
            myfile << Q[i][j] << '\t';
        }
        myfile << '\n';
    }
}

void saveRewardFile() {
	std::ofstream myfile;
	myfile.open("R.txt");
    
	for (int i = 0; i < NBEL; i++) {
		for (int j = 0; j < Nactions; j++) {
			myfile << Rewards[i][j] << '\t';
		}
        myfile	<< '\n';
	}
    
    myfile.close();
    
}

void initialize()
{	int i,j,ii,jj,k;
	int cont;
	int iStdev,iDV;
	//float GaussCumul[NDV][NSTEPS];
	float temp;
	int contador;
    
    //vector of coherences:
    for (i=0;i<NCOH;i++) {
        coherences[i] = metaCoherences[i]*deltatime*KAPPA;
        
        // prior
        p_prior_coh[i] = 1.0;
    }

	// init mu vector
	for (i=0;i<NDV;i++){
		dv[i]=(float) -MAXDV + DVDELTA*i;
	}
	
	
	// init beliefs matrix
	cont = -1;
    for (k=0;k<2;k++) {
        for (i=0;i<NDV;i++){
            for (j=0;j<NSTEPS;j++){
                for (ii=0;ii<NDV;ii++) {
                    for (jj=0;jj<NSTEPS;jj++){
                
                        cont = cont+1;
                        Beliefs[cont][0] = i; //DV1
                        Beliefs[cont][1] = j; //Time step 1
                        Beliefs[cont][2] = ii; //DV2
                        Beliefs[cont][3] = jj; //Time step 2
                        Beliefs[cont][4] = k; //last listened dim
                    }
                }
            }
        }
    }
	
    for (i=0;i<NBEL;i++){
		V[i] 		= 0.0;
		Vprev[i] 	= 0.0;
		
        //		updateOrder[i] = NBEL-1-i;//backwards
        //        updateOrder[i] = i;
	}
    
    //ordeno de mayor iStdev a menor
    cont = -1;
    j = (int) 2*NSTEPS;
    for (;j>=0;j=j-1) {
        //        std::cout << j << std::endl;
        for (i=0;i<NBEL;i++)
        {
            if ((Beliefs[i][1] + Beliefs[i][3])==j)
            {
                cont = cont+1;
                updateOrder[cont] = i;
                
            }
        }
    }
    
}
    


void makeRewardMatrix() {
    int i,j,k,iDV1,iStdev1,iDV2,iStdev2, iLast;
    float pplus1,pplus2,p1,p2,p3,p4, p_diff, w;
    float p_coh1[NCOH],p_coh2[NCOH];
    int method_reward;
    float t_err;
    
    for (i=0;i<NBEL;i++){
		
        //std::cout << i << '\t' << NBEL << '\n';
        
        iDV1 		= Beliefs[i][0];
		iStdev1  	= Beliefs[i][1];
        iDV2        = Beliefs[i][2];
        iStdev2     = Beliefs[i][3];
        iLast       = Beliefs[i][4];
        
        pplus1     = calc_pplus(dv[iDV1], iStdev1, coherences);
        pplus2     = calc_pplus(dv[iDV2], iStdev2, coherences);
        
        
        
        //pminus 	= (float) 1.0-pplus;
        
//        std::cout << pplus << '\t';
        
       
            
        calc_pcoh(p_coh1,dv[iDV1],iStdev1,coherences);
        calc_pcoh(p_coh2,dv[iDV2],iStdev2,coherences);
        
        p_diff = 0.0;
        for (k=0;k<NCOH;k++)
        {
            for (j=0;j<NCOH;j++)
            {
                if (fabs(coherences[k])>fabs(coherences[j]))
                    p_diff += (float) p_coh1[k]*p_coh2[j];
                else if (fabs(coherences[k])==fabs(coherences[j]))
                    p_diff += (float) p_coh1[k]*p_coh2[j]/2.0;
            }
        }
        std::cout << p_diff << '\n';
        
//        w = 0.5; // weight given to the dir choice. If = 0: no weight to it
//
//        p1 = (pplus1)*(p_diff) + (1.0-w)*(1-pplus1)*(p_diff);
//        p2 = (1-pplus1)*(p_diff) + (1.0-w)*(pplus1)*(p_diff);
//        p3 = (pplus2)*(1-p_diff) + (1.0-w)*(1-pplus2)*(1-p_diff);
//        p4 = (1-pplus2)*(1-p_diff) + (1.0-w)*(pplus2)*(1-p_diff);
//
        
        
        method_reward = 1;
        if (method_reward==0) { //usual one
//            Rewards[i][0] = p1*rewardSuccess + (1-p1)*rewardFailure;
//            Rewards[i][1] = p2*rewardSuccess + (1-p2)*rewardFailure;
//            Rewards[i][2] = p3*rewardSuccess + (1-p3)*rewardFailure;
//            Rewards[i][3] = p4*rewardSuccess + (1-p4)*rewardFailure;
//
//            Rewards[i][4] = rewardListen;
//            Rewards[i][5] = rewardListen;
            
        }
        
        
        else if (method_reward==1) { //penalizing time
            //t_err = (1.0-p_diff + 1.0-pplus1 + 1.0-pplus2)*Tpenalty; //not correct
            //t_err = 0.0;
            t_err = (1.0-p_diff)*Tpenalty; //not correct
            
            Rewards[i][0] = p_diff*Rdiff + pplus1*Reasy + pplus2*Rhard  - (TnonDecision+Tintertrial+t_err)*rho;
            Rewards[i][1] = p_diff*Rdiff + pplus1*Reasy + (1.0-pplus2)*Rhard - (TnonDecision+Tintertrial+t_err)*rho;
            Rewards[i][2] = p_diff*Rdiff + (1.0-pplus1)*Reasy + pplus2*Rhard - (TnonDecision+Tintertrial+t_err)*rho;
            Rewards[i][3] = p_diff*Rdiff + (1.0-pplus1)*Reasy + (1.0-pplus2)*Rhard - (TnonDecision+Tintertrial+t_err)*rho;
            
            t_err = (p_diff)*Tpenalty; //not correct
            Rewards[i][4] = (1.0-p_diff)*Rdiff + pplus2*Reasy + pplus1*Rhard - (TnonDecision+Tintertrial+t_err)*rho;
            Rewards[i][5] = (1.0-p_diff)*Rdiff + pplus2*Reasy + (1.0-pplus1)*Rhard - (TnonDecision+Tintertrial+t_err)*rho;
            Rewards[i][6] = (1.0-p_diff)*Rdiff + (1.0-pplus2)*Reasy + pplus1*Rhard - (TnonDecision+Tintertrial+t_err)*rho;
            Rewards[i][7] = (1.0-p_diff)*Rdiff + (1.0-pplus2)*Reasy + (1.0-pplus1)*Rhard - (TnonDecision+Tintertrial+t_err)*rho;
            
            
            if (force_alternation) // and ignore switch costs
            {
                if (iLast==0){//id last sampled
                    Rewards[i][8] = rewardListen - rho*deltatime + DummyPenalty;
                    Rewards[i][9] = rewardListen - rho*deltatime;
                }
                else
                {
                    Rewards[i][8] = rewardListen - rho*deltatime;
                    Rewards[i][9] = rewardListen - rho*deltatime + DummyPenalty;
                }
                
            }
            else 
            {
                if (iLast==0){//id last sampled
                    Rewards[i][8] = rewardListen - rho*deltatime;
                    Rewards[i][9] = rewardListen - rho*(deltatime + TswitchDim);
                }
                else
                {
                    Rewards[i][8] = rewardListen - rho*(deltatime + TswitchDim);
                    Rewards[i][9] = rewardListen - rho*deltatime;
                }
            }
            
            
        }
        
        
        //        std::cout << Rewards[i][0] << '\t';
        //std::cout << Rewards[i][0] << std::endl;
        if (iStdev1==NSTEPS-1)
            Rewards[i][8] = rewardTimesUp - rho*Tpenalty;
        
        if (iStdev2==NSTEPS-1)
            Rewards[i][9] = rewardTimesUp - rho*Tpenalty;
        
        for (j=0;j<8;j++)
            if (Rewards[i][j]!=Rewards[i][j]) { //es nan
                estadoNan[i] = 1;
            }
        else {
            estadoNan[i] = 0;
        }
        
	}
}

void makeTListenVal(){
    //float M[NDV][NSTEPS][NDV];
    float p_coh[NCOH];
    float suma;
    int i,j,k,m;
    for (i=0;i<NDV;i++)
        for (j=0;j<NSTEPS;j++)
        {
            suma = 0.0;
            for (k=0;k<NDV;k++)
            {
                calc_pcoh(p_coh, dv[i], j, coherences);
                TlistenVal[i][j][k] = 0.0;
                for (int ii=0;ii<NCOH;ii++){
                    TlistenVal[i][j][k] +=  (float) p_coh[ii]*normalpdf(dv[k],dv[i]+coherences[ii],SIGMA);
                }
                suma += TlistenVal[i][j][k];
                
            }
            for (k=0;k<NDV;k++)
                TlistenVal[i][j][k] = TlistenVal[i][j][k] / suma;
        }
}


void makeTransitionMatrices(){
	int cont;
	int i,j,ib;
	int iDV1,iStep1,iDV2,iStep2;
	float suma;
	float p_coh[4];
    int index;
	
    makeTListenVal();
    
 	for (i=0;i<NBEL;i++){
		
		iDV1   = Beliefs[i][0];
		iStep1 = Beliefs[i][1];
        iDV2   = Beliefs[i][2];
        iStep2 = Beliefs[i][3];
    
		// Listen Action
        for (j=0;j<NDV;j++){ //j is iDV_next
            
            // if listening to the first dimension
            index = (int)(j)*(NSTEPS*NDV*NSTEPS) + (iStep1+1)*(NDV*NSTEPS) + iDV2*NSTEPS + iStep2;
            TlistenInd[i][j][0] = index;
            
            // if listening to the second dimension
            index = (int) (j)*(NSTEPS) + (iStep2+1) + NSTEPS*NDV*NSTEPS*iDV1 + NDV*NSTEPS*iStep1 + NSTEPS*NDV*NSTEPS*NDV;
            // the last addend is bevause the after this step the last action would be the listening of the 2nd dim
            TlistenInd[i][j][1] = index;
            
        }
        
	}
    
	//start
	cont = 0;
	for (i=0;i<NBEL;i++){
		if (dv[Beliefs[i][0]]==0 & Beliefs[i][1]==0 &
            dv[Beliefs[i][2]]==0 & Beliefs[i][3]==0 &
            Beliefs[i][4]==0){
            TStartInd = i; // Se usa para algo???????
            TStartVal = 1.0; // Se usa para algo???????
			cont = cont+1; //deberia ser 1!!
		}
	}
	
}


double GetMax(float dArray[], int iSize, int *indexOfMax) {
    int iCurrMax = 0;
    float currMax;
    currMax = dArray[0];
    for (int i = 0; i < iSize; ++i) {
        if (dArray[i]>currMax) {
            iCurrMax = i;
            currMax  = dArray[iCurrMax];
        }
    }
    (*indexOfMax) = iCurrMax;
    return currMax;
}




// A utility function to swap two integers
void swap (int *a, int *b)
{
    int temp = *a;
    *a = *b;
    *b = temp;
}



// A function to generate a random permutation of arr[]
void randomize (int arr[], int n)
{
    // Use a different seed value so that we don't get same
    // result each time we run this program
    srand ( time(NULL) );
    
    // Start from the last element and swap one by one. We don't
    // need to run for the first element that's why i > 0
    for (int i = n-1; i > 0; i--)
    {
        // Pick a random index from 0 to i
        int j = rand() % (i+1);
        
        // Swap arr[i] with the element at random index
        swap(&arr[i], &arr[j]);
    }
}


int main()
{
	int i,j,iter;
	int cont;
	int iBelief, iAction;
	double maxVal;
	double delta;
	std::string separador = "-------";
	int ind;
	double p;
    int iDV,iStep;
	int iDV1,iStep1,iDV2,iStep2,dim_to_listen;
//    int iStepCurrent;
//    float dvCurrent, sigmaCurrent;
	int iCurrMax;
	float debugSuma;
    float g;
    
    float rho1,rho2;
	
    rho1 = -40;//supongo asociado con V[init] positivo
    rho2 = 40;//supongo asociado con V[init] negativo
    //	delta = 10;
    V[TStartInd] = maxError+1.0;//para que entre
	while (fabs(V[TStartInd])>maxError) {
        
        rho = (rho1+rho2)/2;
        
        initialize();
        makeRewardMatrix();
        makeTransitionMatrices();
        
        delta = 10;
        while (delta>maxError) {
            //randomize(updateOrder,NBEL);
            
            for (iter=0;iter<NBEL;iter++){
                iBelief = updateOrder[iter];
                
                iDV1   = Beliefs[iBelief][0];
                iStep1 = Beliefs[iBelief][1];
                iDV2   = Beliefs[iBelief][2];
                iStep2 = Beliefs[iBelief][3];
                
                
                //dvCurrent = dv[iDV];
                //sigmaCurrent = sigmaDV[iStep];
                //iStepCurrent = iStep;
                
                for (iAction=0;iAction<Nactions;iAction++){
                    
                    //SetTransitions(&T[0],iBelief,iAction);
                    qact[iAction] = Rewards[iBelief][iAction];
                    
                    //debugSuma = 0;
                    if ((iAction==8 & iStep1<(NSTEPS-1)) || (iAction==9 & iStep2<(NSTEPS-1)))
                    { //listen
                        dim_to_listen = (int) iAction-8; //0 or 1 for first/second dimension
                        
                        if (dim_to_listen==0)
                        {
                            iDV = iDV1;
                            iStep = iStep1;
                        }
                        else
                        {
                            iDV = iDV2;
                            iStep = iStep2;
                        }
                            
                        
                        for (j=0;j<NDV;j++){
                            ind =  TlistenInd[iBelief][j][dim_to_listen];
                            p   = TlistenVal[iDV][iStep][j];
                            //p   =  TlistenVal[iBelief][j];
                            
                            if (estadoNan[ind]==0) {
                                qact[iAction] = qact[iAction] + p*(gama*V[ind]);
                            }
                        }
                    }
                    
                    Q[iBelief][iAction] = qact[iAction]; //only for saving
                    
                }
                // obtengo la mejor accion, y la uso para updatear el valor
                maxVal = GetMax(qact,Nactions,&iCurrMax);
                V[iBelief] = maxVal;
                BestAction[iBelief] = iCurrMax;
                
                //~ for (j=0;j<Nactions;j++){
                //~ std::cout << qact[j] << '\n';
                //~ }
                //~ std::cout << separador << '\n';
                
            }
            
            //for av rew: resto el valor del un estado arbitrario
            //        g = V[TStartInd];
            //		for (i=0;i<NBEL;i++){
            //            V[i] = V[i]-g;
            //        }
            
            //std::cout << V[TStartInd] << std::endl;
            
            delta = 0;
            for (i=0;i<NBEL;i++){
                if (V[i]==V[i]) {
                    Vdelta[i] = fabs(V[i] - Vprev[i]);
                    Vprev[i]  = V[i];
                }
                else { //el valor es nan
                    Vdelta[i] = 0;
                    Vprev[i] = 0;
                }
            }
            
            delta = GetMax(Vdelta,NBEL,&iCurrMax);
            std::cout << delta << '\n';
            
        }
        
        std::cout << V[TStartInd] << '\n';
        
        std::cout << rho << '\n';
//        std::cout << maxError << '\n';
        std::cout << separador << '\n';
        
        
        if (V[TStartInd]>0) {
            rho1 = rho;
        }
        else rho2 = rho;
        
	}
	
	//save results to file
    
    savePars();
	saveToFile();
    // I don't save this for now
    //saveRewardFile();
	//saveQ();
	
	
    //~
	//~
	//~ cont = 0;
    //~ for (i=0;i<NBEL;i++){
    //~ if (T[i]>0) {
    //~ cont=cont+1;
    //~ std::cout << i << '\t' << Beliefs[i][0] << '\t' << Beliefs[i][1] << '\t'<< Beliefs[i][2] << '\t'<< Beliefs[i][3] << '\t'<< Beliefs[i][4] << '\n';
    //~ }
	//~ }
	//~ 
	//~ std::cout << cont << '\n';
	
	return 0;
}
