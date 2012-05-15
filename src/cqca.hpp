#ifndef __CQCA_HPP__
#define __CQCA_HPP__

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <ctime>
#include "qca.hpp"
#include "version.hpp"

// TODO: at least minimal documentation on how Configurator and VConfiguration
//       work

// rtree = result tree
typedef boost::property_tree::basic_ptree<std::string, double> rtree;

class ConfigurationException: public std::runtime_error
{
public:
    ConfigurationException(const std::string& msg)
        : std::runtime_error(msg) {}
};

namespace ptreeHelpers
{
    using boost::property_tree::ptree;

    template<typename T>
    std::vector<T> getArray (const ptree& c)
    {
        std::vector<T> v;
        for (ptree::const_iterator i=c.begin(); i!=c.end(); i++)
            v.push_back(i->second.get_value<T>());
        return v;
    }

    template <class Tree, typename T>
    Tree constructArray (const std::vector<T>& v)
    {
        typedef typename Tree::value_type pair;
        Tree c;
        for (size_t i=0; i<v.size(); i++)
        {
            Tree p;
            p.put_value(v[i]);
            c.push_back(pair("", p));
        }
        return c;
    }

    template<class Tree>
    bool isArray (const Tree& c)
    {
        if (c.data() != typename Tree::data_type())
            return false;
        if (c.size() == 0)
            return false;
        bool flag = true;
        for (typename Tree::const_iterator i=c.begin(); i!=c.end(); i++)
            if (i->first != "") flag = false;
        return flag;
    }

    template<class Tree>
    bool hasKey (const Tree& t, const typename Tree::key_type& k)
    {
        typename Tree::const_assoc_iterator i = t.find(k);
        return i != t.not_found();
    }

    template<typename T>
    T fromString (const std::string& v)
    {
        T returnValue;
        std::stringstream s(v);
        s >> returnValue;
        if (!s.eof() || v.empty())
            throw ConfigurationException("Can't convert '" + v + "'.");
        return returnValue;
    }

    template<typename T>
    std::string toString (const T& v)
    {
        std::stringstream s;
        s << v;
        return s.str();
    }
}

using boost::property_tree::ptree;
using boost::property_tree::ptree_bad_path;
using namespace ptreeHelpers;

class CLayout
{
private:
    enum LayoutType {wire, nonuniformwire};

    LayoutType type;
    int cells;
    double a, b, Pext;
    std::vector<double> bs;
    ElectronsPerCell epc;
    Layout l;
    ptree oc;

public:
    CLayout ()
    {}

    void setConfig (const ptree& c)
    {
        // read in parameters
        type = getType(c);
        cells = c.get("cells", 1);
        a = getA(c);
        b = getB(c);
        bs = getBs(c);
        if (bs.size() == 0) 
            for (int i=0; i<cells; i++) 
                bs.push_back(3.0);
        Pext = c.get("Pext", 0.0);
        epc = getEpc(c);

        // validate parameter values
        if (cells<1)
            throw ConfigurationException("A minimum of 1 cell is required");
        if (Pext<-1 || Pext>1)
            throw ConfigurationException("Pext must be in the range [-1:1]");
        bool bsPos=true;
        for (size_t i=0; i<bs.size(); i++)
            if (bs[i]<0) bsPos=false;
        if (a<0 || b<0 || !bsPos)
            throw ConfigurationException("Wire dimensions (a, b, bs) must be positive");
        if (cells != static_cast<int>(bs.size()))
            throw ConfigurationException("The number of b-values in 'bs' must be"
                        " equal to the number of cells for the non-uniform wire.");

        // construct the wire
        if (type == wire)
            l.wire(cells, a, b, Pext, epc);
        else if (type == nonuniformwire)
            l.nonuniformWire(cells, a, bs, Pext, epc);

        oc = c;
    }

    ptree getConfig () const
    {
        ptree c = oc;

        if (type == wire)
        {
            c.put("type", "wire");
            c.put("cells", cells);
            setAOrV1(c);
            setBOrBoa(c);
            c.put("Pext", Pext);
            c.put("epc", epc);
        }
        else if (type == nonuniformwire)
        {
            c.put("type", "nonuniformwire");
            c.put("cells", cells);
            setAOrV1(c);
            setBsOrBoas(c);
            c.put("Pext", Pext);
            c.put("epc", epc);
        }
        return c;
    }

    Layout& layout ()
    {
        return l;
    }

private:
    LayoutType getType (const ptree& c) const
    {
        std::string typeString = c.get("type", "wire");
        if (typeString == "wire")
            return wire;
        else if(typeString == "nonuniformwire")
            return nonuniformwire;
        else
            throw ConfigurationException("Unknown layout type: '" + typeString + "'");
    }

    ElectronsPerCell getEpc (const ptree& c) const
    {
        int epcInt = c.get("epc", 2);
        if (epcInt == 2)
            return epc2;
        else if(epcInt == 6)
            return epc6;
        else
            throw ConfigurationException("Electrons per Cell (epc) must be either 2 or 6");
    }

    double getA (const ptree& c) const
    {
        if (hasKey(c, "a") && hasKey(c, "V1"))
            throw ConfigurationException("Either specify 'a' or 'V1', but "
                                         "not both at the same time");
        if (hasKey(c, "V1"))
            return 1.0 / c.get<double>("V1");
        return c.get<double>("a", 1.0);
    }
    
    double getB (const ptree& c) const
    {
        if (hasKey(c, "b") && hasKey(c, "boa"))
            throw ConfigurationException("Either specify 'b' or 'boa', but "
                                         "not both at the same time");
        if (hasKey(c, "boa"))
        {
            double a_ = getA(c);
            return c.get<double>("boa") * a_;
        }
        return c.get<double>("b", 3.0);
    }
    
    std::vector<double> getBs (const ptree& c) const
    {
        if (hasKey(c, "bs") && hasKey(c, "boas"))
            throw ConfigurationException("Either specify 'bs' or 'boas', but "
                                         "not both at the same time");
        if (hasKey(c, "boas"))
        {
            double a_ = getA(c);
            std::vector<double> bs_ = getArray<double>(c.get_child("boas"));
            for (size_t i=0; i<bs_.size(); i++)
                bs_[i] *= a_;
            return bs_;
        }
        return getArray<double>(c.get_child("bs", ptree()));
    }

    void setAOrV1 (ptree& c) const
    {
        if (hasKey(c, "V1"))
            c.put("V1", 1.0 / a);
        else
            c.put("a", a);
    }
    
    void setBOrBoa (ptree& c) const
    {
        if (hasKey(c, "boa"))
            c.put("boa", b / a);
        else
            c.put("b", b);
    }

    void setBsOrBoas (ptree& c) const
    {
        if (hasKey(c, "boas"))
        {
            std::vector<double> boas = bs;
            for (size_t i=0; i<boas.size(); i++)
                boas[i] /= a;
            c.put_child("boas", constructArray<ptree>(boas));
        }
        else
            c.put_child("bs", constructArray<ptree>(bs));
    }
};

class CQcaGenericBase
{
public:
    virtual ~CQcaGenericBase () {}
    virtual void setConfig (const ptree& c) = 0;
    virtual ptree getConfig () const = 0;
    virtual rtree measure () = 0;
};

template <class QcaSystem>
class CQcaGeneric : public CQcaGenericBase
{
private:
    CLayout l;
    ptree os;
    QcaSystem s;
    ptree oc;

public:
    virtual ~CQcaGeneric () {}

    virtual void setConfig (const ptree& c)
    {
        s.t = c.get("t", 1.0);
        s.td = c.get("td", 0.0); 
        s.Vext = c.get("Vext", 0.0);
        s.V0 = c.get("V0", 1000.0); 
        s.mu = c.get("mu", 0.0);
        s.epsilonr = c.get("epsilonr", QCA_NATURAL_EPSILON_R);
        s.lambdaD = c.get("lambdaD", 0.0);
        s.q = c.get("q", 0);
        s.beta = getBeta(c);
        os = c.get_child("observables", ptree());
        l.setConfig(c.get_child("layout", ptree()));
        s.l = l.layout();

        oc = c;
    }

    virtual ptree getConfig () const
    {
        /*
         * Ideally we would read the layout configuration back from s.l. However, the
         * the implementation of the class Layout currently does not know about its
         * high-level layout (e.g. is it a wire or a nonuniform wire). That's
         * why we store an instance of CLayout in CQca and use CLayout to get
         * the layout configuration.
         */
        ptree c = oc;
        c.put("t", s.t);
        c.put("td", s.td);
        c.put("Vext", s.Vext);
        c.put("V0", s.V0);
        c.put("mu", s.mu);
        c.put("epsilonr", s.epsilonr);
        c.put("lambdaD", s.lambdaD);
        c.put("q", s.q);
        setBetaOrT(c, s.beta);
        c.put_child("layout", l.getConfig());
        c.put_child("observables", os);
        return c;
    }

    virtual rtree measure ()
    {
        //TODO: only call update when necessary -- e.g. not necessary for
        //changed beta
        s.update();

        rtree r;
        for (ptree::const_iterator i=os.begin(); i!=os.end(); i++)
        {
            if (i->first == "P")
                r.put_child("P", measureP(i->second));
            else if (i->first == "P2")
                r.put_child("P2", measureP2(i->second));
            else if (i->first == "N")
                r.put_child("N", measureN(i->second));
            else if (i->first == "E")
            {
                if (i->second.get_value<std::string>() == "yes")
                    r.put_child("E", measureE(i->second));
                else if (i->second.get_value<std::string>() != "no")
                    throw ConfigurationException("For observable 'E': "
                                    "value must be either 'yes' or 'no'");
            }
            else
                throw ConfigurationException("Unknown observable: '" + i->first + "'");
        }
        return r;
    }

    QcaSystem& system ()
    {
        return s;
    }

private:
    std::vector<int> getCells (const ptree& c) const
    {
        std::vector<int> allCells;
        for (int i=0; i<s.N_p(); i++)
            allCells.push_back(i);

        std::vector<int> cells;
        if (c.get_value<std::string>() == "all")
            cells = allCells;
        else if (isArray(c))
            cells = getArray<int>(c);
        else if (c.get_value<std::string>() != "")
            try
            {
                cells.push_back(c.get_value<int>());
            }
            catch (std::runtime_error)
            {
                throw ConfigurationException("Invalid cell specification for observable");
            }

        for (size_t i=0; i<cells.size(); i++)
        {
            if (cells[i] < 0)
                throw ConfigurationException("In observable specification: "
                                             "Cell numbers must be positive");
            if (cells[i] >= s.N_p())
                throw ConfigurationException(
                        "In observable specification: Trying to measure cell " + 
                        toString(cells[i]) + ", but there are only " + 
                        toString(s.N_p()) + " cells in the system.");
        }

        return cells;
    }

    rtree measureP (const ptree& c) const
    {
        std::vector<int> cells = getCells(c);
        rtree r;
        for (size_t i=0; i<cells.size(); i++)
            r.put(toString(cells[i]), s.measurePolarization(cells[i]));
        return r;
    }

    rtree measureP2 (const ptree& c) const
    {
        std::vector<int> cells = getCells(c);
        rtree r;
        for (size_t i=0; i<cells.size(); i++)
            r.put(toString(cells[i]), s.measurePolarization2(cells[i]));
        return r;
    }

    rtree measureN (const ptree& c) const
    {
        std::vector<int> cells = getCells(c);
        rtree r;
        for (size_t i=0; i<cells.size(); i++)
        {
            std::vector<double> ns = s.measureParticleNumber(cells[i]);
            r.put(toString(cells[i]) + ".total", ns[4]);
            r.put(toString(cells[i]) + ".0", ns[0]);
            r.put(toString(cells[i]) + ".1", ns[1]);
            r.put(toString(cells[i]) + ".2", ns[2]);
            r.put(toString(cells[i]) + ".3", ns[3]);
        }
        return r;
    }

    class EpsilonLess
    {
    public:
        EpsilonLess (double e_ = 1E-20)
        : e(e_)
        {}

        bool operator() (double a, double b)
        {
            return b-a > e;
        }
    private:
        double e;
    };

    rtree measureE (const ptree& c) 
    {
        rtree r;
        typedef std::map<double, int, EpsilonLess> myMap;
        typedef typename myMap::const_iterator mapIt;
        myMap evs;
        for (int i=0; i<s.energies().size(); i++)
            evs[s.energies()(i)]++;
        int j=0;
        for (mapIt i=evs.begin(); i!=evs.end(); i++)
        {
            rtree v;
            v.put("abs", i->first);
            v.put("rel", i->first - s.Emin());
            v.put("deg", static_cast<double>(i->second));
            r.push_back(rtree::value_type(toString(j), v));
            j++;
        }
        return r;
    }

    double getBeta (const ptree& c) const
    {
        if (hasKey(c, "beta") && hasKey(c, "T"))
            throw ConfigurationException("Either specify 'beta' or 'T', but "
                                         "not both at the same time");
        if (hasKey(c, "T"))
            return 1.0 / c.get<double>("T");
        return c.get<double>("beta", 1);
    }

    void setBetaOrT (ptree& c, double beta) const
    {
        if (hasKey(c, "T"))
            c.put("T", 1.0 / beta);
        else
            c.put("beta", beta);
    }
};

class CQca
{
private:
    enum ModelType {none, bond, fixed, grand};
    
    CQcaGenericBase* s;
    ModelType model;

public:
    CQca ()
    : model(none)
    {}

    ~CQca ()
    {
        if (model != none)
            delete s;
    }
    
    void setConfig (const ptree& c)
    {
        ModelType nm = getModel(c.get("model", "grand"));
        if (model != nm)
        {
            if (model != none)
                delete s;
            model = nm;
            if (model == bond)
                s = new CQcaGeneric<QcaBond>();
            else if (model == fixed)
                s = new CQcaGeneric<QcaFixedCharge>();
            else if (model == grand)
                s = new CQcaGeneric<QcaGrandCanonical>();
        }
        s->setConfig(c);
    }

    ptree getConfig () const
    {
        if (model == none)
            return ptree();
        // it's nice to have "model" as the first entry in the configuration,
        // that's why we're using c.push_front instead of the simpler c.put()
        ptree c = s->getConfig();
        if (c.find("model") == c.not_found())
            c.push_front(ptree::value_type("model", ptree(modelToString(model))));
        return c;
    }

   rtree measure ()
   {
       if (model == none)
           return rtree();
       return s->measure();
   }

private:
   ModelType getModel (const std::string& s) const
   {
       if (s == "bond")
           return bond;
       if (s == "grand" || s == "grandcanonical")
           return grand;
       if (s == "fixed" || s == "fixedcharge")
           return fixed;
       throw ConfigurationException("Invalid model: '" + s + "'");
   }

   std::string modelToString (const ModelType m) const
   {
       if (m == bond)
           return "bond";
       if (m == fixed)
           return "fixedcharge";
       if (m == grand)
           return "grandcanonical";
       return "none";
   }
};

class VParam
{
/*
 * This implementation is not very good. 
 * 1. It really is two classes with the same interface in one (generator and
 *    normal vparam).
 * 2. We rely on a pointer to the "global" configuration ptree c. So
 *    VParam only works as long as this configuration variable remains valid.
 */
private:
    std::string name_;
    bool IAmAnArray;
    
    // Has to be a pointer. Otherwise it will be invalid when this instance of
    // VParam is copied (happens, for example, if we use VParam in a vector)
    ptree* t;
    ptree::const_iterator it;

    std::vector<double> values;
    int index;

public:
    VParam (ptree&  c, const std::string& name__)
    : name_(name__)
    {
        try {
            c.get_child(name_);
        }
        catch (ptree_bad_path e)
        {
            throw ConfigurationException("Changing parameter '" + name_ + 
                                         "' does not exist.");
        }

        ptree& p = c.get_child(name_);
        if (isArray(p))
        {
            IAmAnArray = true;
            t = &p;
            it = t->begin();
        }
        else if (isGenerator(p))
        {
            IAmAnArray = false;
            values = getArrayFromGenerator(p);
            index = 0;
        }
        else
            throw ConfigurationException("Changing parameter '" + name_ + 
                                         "' has invalid specification.");
    }

    ptree value () const
    {
        if (IAmAnArray)
            return it->second;
        else
        {
            ptree n;
            n.put_value(toString(values[index]));
            return n;
        }
    }

    size_t size () const
    {
        if (IAmAnArray)
            return t->size();
        else
            return values.size();
    }

    void increment () 
    {
        if (IAmAnArray)
        {
            it++;
            if (it == t->end())
                it = t->begin();
        }
        else
        {
            index++;
            if (index == values.size())
                index = 0;
        }
    }

    bool atFirstElement () const
    {
        if (IAmAnArray)
            return it == t->begin();
        else
            return index == 0;
    }

    std::string name () const
    {
        return name_;
    }

private:
    bool isGenerator (const ptree& c) const
    {
        return (c.size() == 0) && (c.data().find(';') != c.data().npos);
    }

    std::vector<double> getArrayFromGenerator (const ptree& c) const
    {
        std::vector<double> tuple = getList<double>(c.data(), ';');
        if (tuple.size() != 3) 
            throw ConfigurationException("Invalid generator specification: '" + c.data() + "'.");
        double begin = tuple[0];
        double end = tuple[1];
        double inc = tuple[2];
        double epsilon = inc * 0.0001;
        if (inc==0 || (begin<=end && inc<0) || (end<=begin && inc>0)) 
            throw ConfigurationException("Invalid generator specification: '" + c.data() + "'.");
        std::vector<double> generated;
        for (double i=begin; (begin<=end && i<=end+epsilon) || 
                             (end<=begin && i>=end+epsilon); i+=inc)
            generated.push_back(i);
        return generated;
    }

    template<typename T>
    std::vector<T> getList (const std::string& v, const char sep) const
    {
        if (v.size() == 0) return std::vector<T>();
        std::vector<T> list;
        size_t pos1=0;
        size_t pos2=0;
        while (pos1 <= v.size())
        {
            pos2 = v.find(sep, pos1);
            if (pos2 == v.npos) pos2 = v.size();
            if (pos2 == pos1) 
                throw ConfigurationException("Error parsing list: '" + v + "'.");
            T item = fromString<T>(v.substr(pos1, pos2-pos1));
            list.push_back(item);
            pos1 = pos2+1;
        }
        return list;
    }
};

class VConfiguration
{
private:
    ptree c;
    bool gotConfig, gotAllVariants, constructedVariants, hasVariants;
    std::vector<VParam> vparams;
    std::vector<std::string> changed;

public:
    VConfiguration (const ptree& c_)
    : c(c_), gotConfig(false), gotAllVariants(false), 
      constructedVariants(false), hasVariants(false)
    {}

    ptree getNext () 
    {
        constructVariants();

        if (!hasNext())
            throw ConfigurationException("This configuration has no more variants.");
        
        if (!hasVariants) 
        {
            gotConfig = true;
            return c;
        }

        ptree p = getConfigForCurrentVariant();
        increaseVariantIndex();
        return p;
    }

    bool hasNext () 
    {
        constructVariants();

        if (!hasVariants) 
            return !gotConfig;

        return !gotAllVariants;
    }

    int numberOfVariants ()
    {
        constructVariants();
        int n=1;
        for (size_t i=0; i<vparams.size(); i++)
            n *= vparams[i].size();
        return n;
    }

private:
    void constructVariants ()
    {
        if (constructedVariants)
            return;
        
        const ptree p = c.get_child("changing", ptree());
        if (!isArray(p))
        {
            if (p.get_value<std::string>() != "")
                throw ConfigurationException("'changing' must be an array of changing parameters");
            hasVariants = false;
            return;
        }

        hasVariants = true;
        changed.clear();
        std::vector<std::string> names = getArray<std::string>(p);
        for (size_t i=0; i<names.size(); i++)
        {
            vparams.push_back(VParam(c, names[i]));
            // initially all varying parameters are considered changed
            changed.push_back(names[i]);
        }
        constructedVariants = true;
    }

    ptree getConfigForCurrentVariant () const
    {
        ptree nc = c;
        for (size_t i=0; i<vparams.size(); i++)
        {
            // replace the changing parameter by it's current value
            // and attach the original config of the changing parameter as
            // "original.parameter_name"
            ptree o = nc.get_child(vparams[i].name());
            nc.put_child("original." + vparams[i].name(), o);
            nc.put_child(vparams[i].name(), vparams[i].value());
        }
        nc.put_child("changed", constructArray<ptree>(changed));
        return nc;
    }

    void increaseVariantIndex ()
    {
        changed.clear();
        for (int i=vparams.size()-1; i>=0; i--)
        {
            VParam& v = vparams[i];
            changed.push_back(v.name());
            // note: increment() wraps around
            // so atFirstElement() checks if we just wrapped around
            v.increment(); 
            if (!v.atFirstElement())
                break;
        }
        // keep original order of entries in changing
        std::reverse(changed.begin(), changed.end());
        
        gotAllVariants = true;
        for (size_t i=0; i<vparams.size(); i++)
            if (!vparams[i].atFirstElement())
            {
                gotAllVariants = false;
                return;
            }
    }
};

class Configurator 
{
private:
    std::vector<VConfiguration> cs;
    size_t i_c;

public:
    Configurator (const std::string& json)
    : i_c(0)
    {
        ptree c;
        std::string js;
        if (isfile(json))
            // could be implemented in a better, more efficient way
            js = getFileAsString(json);
        else
            js = json;
        std::stringstream ss(jsonify(js));
        read_json(ss, c);
        
        if (isArray(c))
            for (ptree::const_iterator i=c.begin(); i!=c.end(); i++)
                cs.push_back(VConfiguration(i->second));
        else
            cs.push_back(VConfiguration(c));
    }

    bool hasNext ()
    {
        if (i_c == cs.size()-1)
            return cs[i_c].hasNext();
        return i_c < cs.size();
    }

    ptree getNext ()
    {
        if (cs[i_c].hasNext())
            return cs[i_c].getNext();
        
        i_c++;
        if (i_c < cs.size())
            return cs[i_c].getNext();
        else
            throw ConfigurationException("No more configurations available.");
    }

    int numberOfConfigs ()
    {
        int sum=0;
        for (size_t i=0; i<cs.size(); i++)
            sum += cs[i].numberOfVariants();
        return sum;
    }

    static std::string jsonify (const std::string& s)
    {
        std::stringstream is(s);
        std::stringstream os;

        bool needClosingBracket = false;
        if (is.peek() != '{' && is.peek() != '[')
        {
            os << "{";
            needClosingBracket = true;
        }

        bool meQuoting = false;
        bool userQuoting = false;
        char c;
        while (is >> std::noskipws >> c)
        {
            if (c == '\'' || c == '"')
            {
                if (userQuoting) userQuoting = false;
                else userQuoting = true;
                os << "\"";
            }
            else if (userQuoting)
                os << c;
            else if (std::isspace(c))
                continue;
            else if (c == '{' || c == '}' || c == '[' || c == ']' || c == ',' || c == ':')
            {
                if (meQuoting)
                {
                    os << "\"";
                    meQuoting = false;
                }
                os << c;
            }
            else if (!meQuoting)
            {
                os << "\"" << c;
                meQuoting = true;
            }
            else
                os << c;
        }

        if (meQuoting)
            os << "\"";
        if (needClosingBracket)
            os << "}";

        return os.str();
    }

    static std::string getFileAsString (const std::string& file)
    {
        std::ifstream is(file.c_str());
        std::stringstream os;
        while (is.good())
        {
            const int size = 200;
            char cs[size];
            is.read(cs, size);
            os.write(cs, is.gcount());
        }
        return os.str();
    }

private:
    bool isfile(std::string path)
    {
        struct stat s;
        if (stat(path.c_str(), &s) != 0) {
            return false;
        }
        if (S_ISREG(s.st_mode)) {
            return true;
        }
        return false;
    }
};

class Store
{
private:
    ptree c_old;
    int lineCount;
    // column width for output; how often (every x lines) to print table header
    enum {columnWidth=14, tableHeaderEveryXLines=20};

    template<class Tree>
    class TreeTypes
    {
    public:
        typedef std::pair<typename Tree::key_type, typename Tree::data_type> pair;
        typedef std::vector<pair> vector;
    };

public:
    Store ()
    : lineCount(0)
    {}

    void store (const ptree& c, const rtree& r)
    {
        // configure output
        std::cout << std::setprecision(6);
        std::cout << std::left;

        ptree c_new = removeVParams(c);
        if (c_new != c_old)
        {
            printHeader(c);
            c_old = c_new;
        }
        printTableRow(c,r);
    }

private:
    ptree removeVParams (ptree c) const
    {
        // replace vparams with an empty node
        std::vector<std::string> names = getVParams(c);
        for (size_t i=0; i<names.size(); i++)
            c.put_child(names[i], ptree());
        c.put_child("changed", ptree());
        return c;
    }

    void printHeader (const ptree& c)
    {
        // only false if this is the first header / measurement
        if (lineCount != 0)
            std::cout << std::endl;

        std::cout << "# program version: " << GIT_PROGRAM_VERSION << std::endl 
                  << "# date: " << getDate() << std::endl
                  << "# " << std::endl;
        
        std::stringstream ss;
        write_json(ss, prepareForPrinting(c));
        std::string s = ss.str();
        prependLines(s, "# ");
        std::cout << s;
        std::cout << "# " << std::endl;
        
        lineCount = 0;
    }

    void printTableHeader (const ptree& c, const rtree& r) const
    {
        std::vector<std::string> ps = getVParams(c);
        TreeTypes<rtree>::vector os = getFlattenedTree(r);
        std::cout << "# ";
        for (size_t i=0; i<ps.size(); i++)
        {
            const ptree& n = c.get_child(ps[i]);
            if (isArray(n))
                for (size_t j=0; j<n.size(); j++)
                    std::cout << std::setw(columnWidth) << (ps[i] + "." + toString(j));
            else 
                std::cout << std::setw(columnWidth) << ps[i];
        }
        for (size_t i=0; i<os.size(); i++)
            std::cout << std::setw(columnWidth) << os[i].first;
        std::cout << std::endl;
    }

    void printTableRow (const ptree& c, const rtree& r)
    {
        if (r.size() == 0)
            return;

        if (lineCount % tableHeaderEveryXLines == 0)
            printTableHeader(c,r);

        std::vector<std::string> ps = getVParams(c);
        TreeTypes<rtree>::vector os = getFlattenedTree(r);
        std::cout << "  ";
        for (size_t i=0; i<ps.size(); i++)
        {
            const ptree& n = c.get_child(ps[i]);
            if (isArray(n))
                for (ptree::const_iterator j=n.begin(); j!=n.end(); j++)
                    std::cout << std::setw(columnWidth) << j->second.get_value<std::string>();
            else
                std::cout << std::setw(columnWidth) << c.get<double>(ps[i], 0);
        }
        for (size_t i=0; i<os.size(); i++)
            std::cout << std::setw(columnWidth) << os[i].second;
        std::cout << std::endl;
        lineCount++;
    }

    std::vector<std::string> getVParams (const ptree& c) const
    {
        ptree p = c.get_child("changing", ptree());
        if (!isArray(p))
            return std::vector<std::string>();
        return getArray<std::string>(p);
    }

    template<class Tree>
    typename TreeTypes<Tree>::vector getFlattenedTree (const Tree& t) const
    {
        typename TreeTypes<Tree>::vector v;
        if (t.size() != 0)
            getFlattenedTreeRecursive(t, v, "");
        return v;
    }
    
    template<class Tree>
    void getFlattenedTreeRecursive (const Tree& t, 
                                    typename TreeTypes<Tree>::vector& v, 
                                    const std::string& path) const
    {
        if (t.size() == 0)
        {
            v.push_back(typename TreeTypes<Tree>::pair(path, t.data()));
            return;
        }

        bool isarray = isArray(t);
        int l=0;
        for (typename Tree::const_iterator i=t.begin(); i!=t.end(); i++)
        {
            std::string npath = path;
            if (npath.length() != 0)
                npath += ".";
            if (isarray)
                npath += toString(l);
            else
                npath += i->first;
            getFlattenedTreeRecursive(i->second, v, npath);
            l++;
        }
    }

    std::string getDate () const
    {
        time_t t;
        tm localTime;
        char cstr[100];
        time(&t);
        tzset();
        localtime_r(&t, &localTime);
        strftime(cstr, 100, "%a %b %d %T %Y %z", &localTime);
        return cstr;
    }
    
    void prependLines (std::string& s, const std::string& pr) const
    {
        size_t pos = 0;
        while(pos<s.length() && pos!=s.npos)
        {
            s.insert(pos, pr);
            pos = s.find('\n', pos);
            pos++;
        }
    }
    
    ptree prepareForPrinting (ptree c)
    {
        ptree::assoc_iterator i;
        i = c.find("changed");
        if (i != c.not_found())
            c.erase(c.to_iterator(i));
        i = c.find("original");
        if (i != c.not_found())
        {
            // replaces the current value of a varying parameter by its
            // original configuration / specification
            std::vector<std::string> ps = getVParams(c);
            for (size_t j=0; j<ps.size(); j++)
                c.put_child(ps[j], c.get_child("original." + ps[j]));
            c.erase(c.to_iterator(i));
        }
        
        return c;
    }
};

class Runner
{
public:
    static void run (const std::string& jsonConfig)
    {
        Configurator c(jsonConfig);
        CQca s;
        Store o;

        while (c.hasNext())
        {
            s.setConfig(c.getNext());
            o.store(s.getConfig(), s.measure());
        }
    }
};

#endif // __CQCA_HPP__
