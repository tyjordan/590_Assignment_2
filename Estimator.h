#ifndef _ESTIMATOR_HEADER_
#define _ESTIMATOR_HEADER_

#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <cassert>
#include <typeinfo>

#include "Particle.h"
#include "Material.h"
#include "Reaction.h"

class estimator {
  private:
    std::string estimator_name;
  protected:
	double m = 1;
    unsigned long long nhist;
  public:
     estimator( std::string label ) : estimator_name(label) {};
    ~estimator() {};

    virtual std::string name() final { return estimator_name; };

	virtual double multiplier() final { return m; };
	virtual void set_multiplier(double M) final { m = M; return; };

    virtual void score( particle*, double ) = 0;

    virtual void endHistory()       = 0;
    virtual void report()           = 0;
};

class single_valued_estimator : public estimator {
  private:

  protected:
    double tally_hist, tally_sum, tally_squared;
  public:

     single_valued_estimator(std::string label ) : estimator(label) { 
       nhist         = 0;
       tally_hist    = 0.0;   
       tally_sum     = 0.0; 
       tally_squared = 0.0;
     };
    ~single_valued_estimator() {};

     virtual void endHistory()    final { 
       nhist++;
       tally_sum     += tally_hist;
       tally_squared += tally_hist * tally_hist;
       tally_hist = 0.0; }

     virtual void score( particle*, double ) = 0;

     virtual void report() final {
       double mean = tally_sum / nhist;
       double var  = ( tally_squared / nhist - mean*mean ) / nhist;
       std::cout << name() << "   " << mean << "   " << std::sqrt( var ) / mean << std::endl;  
     };

};

class surface_current_estimator : public single_valued_estimator {
  private:

  public:
     surface_current_estimator( std::string label ) : single_valued_estimator(label) {};
    ~surface_current_estimator() {};

    void score( particle*, double );
};

class cell_pathLengthFlux_estimator : public single_valued_estimator {
  public:
     cell_pathLengthFlux_estimator( std::string label ) : 
       single_valued_estimator(label) {};
    ~cell_pathLengthFlux_estimator() {};

    void score( particle*, double );
};

class counting_estimator : public estimator {
  private:
    int count_hist;
    std::vector< double > tally;
  public:
     counting_estimator( std::string label ) : estimator(label) { count_hist = 0; };
    ~counting_estimator() {};

    void score( particle*, double);
    void endHistory();
    void report();
};

#endif
