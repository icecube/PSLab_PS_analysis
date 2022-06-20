#ifndef LLH_MINUITWRAPPERCLASS_H_
#define LLH_MINUITWRAPPERCLASS_H_


// TMinuit has to point to a function to minimize.  Unfortunately, this 
// function cannot be a member function of a class, which would have made it
// easy to access all of the data needed for the minimization.

// The solution is to have a generic function which uses 
// a static external pointer, and to have the static external pointer 
// set to point to your particular data object.

// You make your data class inherit from the MinuitWrapperClass, 
// and define the member function EvalMinuitFCN  to do exactly what you want
// with your data.



class MinuitWrapperClass {
public:
  virtual ~MinuitWrapperClass() { }

  // the base version does nothing!  Must define for you own derived class,
  // if you plan on using minuit
  virtual void EvalMinuitFCN(int &npar, double *gin, double &f, 
			     double *par, int iflag)
  {
    // touch these variables just to eliminate warnings during compile    
    if (npar || gin || f || par || iflag) { }
  };

  void SetMinuitWrapperPtr(MinuitWrapperClass* ptr);
};



static MinuitWrapperClass* ptrMinuitWrapperClass_ = NULL;



void MinuitWrapperFCN (int &npar, double *gin, double &f, 
				 double *par, int iflag)
{
  assert(ptrMinuitWrapperClass_);
  ptrMinuitWrapperClass_->EvalMinuitFCN(npar, gin, f, par, iflag);
}



void MinuitWrapperClass::SetMinuitWrapperPtr(MinuitWrapperClass* ptr) {
  ptrMinuitWrapperClass_ = ptr;
}



/*    

//   EXAMPLE:


// MyClass INHERITS FROM MinuitWrapperClass

class MyClass : MinuitWrapperClass {  
private:
  double myData;

public:

  // THIS IS THE FUNCTION YOU _WISH_ YOU COULD PASS TO MINUIT

  virtual void EvalMinuitFCN(int &npar, double *gin, double &f, 
			     double *par, int iflag) 
  {
    f = myData*myData;
    // etc.   You can use class member data to evaluate f
  }


  // THIS IS HOW YOUR MINIMIZATION ROUTINE IS SET UP

  void DoMinimiziation() {
    TMinuit m;
    // initialize m parameters, then

    m.SetFCN(MinuitWrapperFCN);
    SetMinuitWrapperPtr(this);      // THE WRAPPER FCN WILL USE _YOUR_ FCN
    m.Migrad();
    SetMinuitWrapperPtr(NULL);      // RESET POINTER
 
}

*/


#endif // LLH_MINUITWRAPPERCLASS_H_
