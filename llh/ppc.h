namespace xppc{
  float xrnd();
  float grnd();

  void start();
  void stop();
  void choose(int);
  void ini();
  void fin();
  void eout();

  void flone(unsigned long long);

  struct mcid:std::pair<int,unsigned long long>{
    int frame;
    double t0;

    bool operator< (const mcid & rhs) const {
      return
	frame!=rhs.frame ? frame<rhs.frame :
	first!=rhs.first ? first<rhs.first :
	second!=rhs.second ? second<rhs.second :
	t0<rhs.t0;
    }
  };

  struct DOM{
    float R, F;
    float r[3];
  };

  struct ikey{
    int str, dom;

    bool isinice() const{
      return str>86 || (str>0 && dom>=1 && dom<=60);
    }

    bool operator< (const ikey & rhs) const {
      return str == rhs.str ? dom < rhs.dom : str < rhs.str;
    }

    bool operator!= (const ikey & rhs) const {
      return str != rhs.str || dom != rhs.dom;
    }
  };

  struct OM:DOM,ikey{};
  extern std::vector<OM> i3oms;

  extern std::map<ikey, int> omts;
  extern std::map<ikey, float> hvs;
  extern std::map<ikey, std::pair<float, int> > rdes;

  void initialize(float);
  const DOM& flset(int, int);
  const float * fldir();
  void flshift(float [], float [], float * = NULL);

  struct ihit{
    mcid track;
    ikey omkey;
    int pmt;
    float dir;
    float time;

    bool operator< (const ihit & rhs) const {
      return
	track!=rhs.track ? track<rhs.track :
	omkey!=rhs.omkey ? omkey<rhs.omkey :
	pmt!=rhs.pmt ? pmt<rhs.pmt :
	dir!=rhs.dir ? dir<rhs.dir :
	time<rhs.time;
    }
  };

  void set_res(float);
  void set_res(float, float);

  struct pout{
    float r[4], n[4];
  };

  typedef std::map<ihit, std::vector<pout> > outz;
  extern outz hitz;
  extern bool nextgen;

  void efin();

  void sett(float, float, float, std::pair<int,unsigned long long>, int, double = 0);
  void addp(float, float, float, float, float, float, int, float = 1);

  void addp_clst(float, float, float, float, unsigned long long, float, float);
  void addp_mopo(float, float, float, float, float, float, float, float, float);
}
