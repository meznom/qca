#ifndef __CQCA_HPP__
#define __CQCA_HPP__

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <ctime>
#include "qca.hpp"
#include "version.hpp"

// TODO: in store output: get rid of 'changed'; output original configuration in
//       header (with ranges; e.g. "b": [1,2,3,4])
// TODO: implement start:end:interval syntax
// TODO: at least minimal documentation on how Configurator and VConfiguration
//       work
// TODO: maybe use modified json -> get rid of '"', and treat all literals as
//       strings => easier to write and read


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
        for (int i=0; i<v.size(); i++)
        {
            Tree p;
            p.put_value(v[i]);
            c.push_back(pair("", p));
        }
        return c;
    }

    bool isArray (const ptree& c)
    {
        if (c.get_value<std::string>() != "")
            return false;
        if (c.size() == 0)
            return false;
        bool flag = true;
        for (ptree::const_iterator i=c.begin(); i!=c.end(); i++)
            if (i->first != "") flag = false;
        return flag;
    }
}

using boost::property_tree::ptree;
using namespace ptreeHelpers;
using namespace boost::property_tree::json_parser;

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

public:
    CLayout ()
    {}

    void setConfig (const ptree& c)
    {
        // read in parameters
        type = getType(c);
        cells = c.get("cells", 1);
        a = c.get("a", 1.0);
        b = c.get("b", 3.0);
        bs = getArray<double>(c.get_child("bs", ptree()));
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
        for (int i=0; i<bs.size(); i++)
            if (bs[i]<0) bsPos=false;
        if (a<0 || b<0 || !bsPos)
            throw ConfigurationException("Wire dimensions (a, b, bs) must be positive");
        if (cells != bs.size())
            throw ConfigurationException("The number of b-values in 'bs' must be"
                        " equal to the number of cells for the non-uniform wire.");

        // construct the wire
        if (type == wire)
            l.wire(cells, a, b, Pext, epc);
        else if (type == nonuniformwire)
            l.nonuniformWire(cells, a, bs, Pext, epc);
    }

    ptree getConfig () const
    {
        ptree c;

        if (type == wire)
        {
            c.put("type", "wire");
            c.put("cells", cells);
            c.put("a", a);
            c.put("b", b);
            c.put("Pext", Pext);
            c.put("epc", epc);
        }
        else if (type == nonuniformwire)
        {
            c.put("type", "nonuniformwire");
            c.put("cells", cells);
            c.put("a", a);
            c.put_child("bs", constructArray<ptree>(bs));
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
        s.beta = c.get("beta", 1);
        os = c.get_child("observables", ptree());
        l.setConfig(c.get_child("layout", ptree()));
        s.l = l.layout();
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
        ptree c;
        c.put("t", s.t);
        c.put("td", s.td);
        c.put("Vext", s.Vext);
        c.put("V0", s.V0);
        c.put("mu", s.mu);
        c.put("epsilonr", s.epsilonr);
        c.put("lambdaD", s.lambdaD);
        c.put("q", s.q);
        c.put("beta", s.beta);
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
            else if (i->first == "E" && i->second.get_value<std::string>() == "yes")
                r.put_child("E", measureE(i->second));
            else
                throw ConfigurationException("Unknown observable: '" + i->first + "'");
        }
        return r;
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

        for (int i=0; i<cells.size(); i++)
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
        for (int i=0; i<cells.size(); i++)
            r.put(toString(cells[i]), s.measurePolarization(cells[i]));
        return r;
    }

    rtree measureP2 (const ptree& c) const
    {
        std::vector<int> cells = getCells(c);
        rtree r;
        for (int i=0; i<cells.size(); i++)
            r.put(toString(cells[i]), s.measurePolarization2(cells[i]));
        return r;
    }

    rtree measureN (const ptree& c) const
    {
        std::vector<int> cells = getCells(c);
        rtree r;
        for (int i=0; i<cells.size(); i++)
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
        for (mapIt i=evs.begin(); i!=evs.end(); i++)
        {
            std::vector<double> v(3);
            v[0] = i->first;
            v[1] = i->first - s.Emin();
            v[2] = i->second;
            r.push_back(rtree::value_type("", constructArray<rtree>(v)));
        }
        return r;
    }

    template<typename T>
    std::string toString (T v) const
    {
        std::stringstream ss;
        ss << v;
        return ss.str();
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

class VConfiguration
{
private:
    struct VParam
    {
        VParam (const std::string& name_, 
                const std::vector<std::string> values_, 
                int index_)
        : name(name_), values(values_), index(index_) 
        {}

        std::string currentValue () const
        {
            return values[index];
        }
        
        std::string name;
        std::vector<std::string> values;
        int index;
    };

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
        for (int i=0; i<vparams.size(); i++)
            n *= vparams[i].values.size();
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
        for (int i=0; i<names.size(); i++)
        {
            const ptree a = c.get_child(names[i], ptree());
            if (!isArray(a))
                throw ConfigurationException("Changing parameter '" + names[i] + 
                                             "' does not exist or is not an array");
            VParam v(names[i], getArray<std::string>(a), 0);
            vparams.push_back(v);
            
            // initially all varying parameters are considered changed
            changed.push_back(names[i]);
        }
        constructedVariants = true;
    }

    ptree getConfigForCurrentVariant () const
    {
        ptree nc = c;
        for (int i=0; i<vparams.size(); i++)
        {
            // first delete all children (the array elements) and then set the
            // single new value
            ptree& n = nc.get_child(vparams[i].name);
            for (ptree::iterator i=n.begin(); i!=n.end(); i = n.erase(i));
            n.put_value(vparams[i].currentValue());
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
            v.index++;
            changed.push_back(v.name);
            if (v.index == v.values.size())
                v.index = 0;
            else
                break;
        }
        // keep original order of entries in changing
        std::reverse(changed.begin(), changed.end());
        
        gotAllVariants = true;
        for (int i=0; i<vparams.size(); i++)
            if (vparams[i].index != 0)
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
    int i_c;

public:
    Configurator (const std::string& jsonString)
    : i_c(0)
    {
        ptree c;
        std::stringstream ss(jsonString);
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
        for (int i=0; i<cs.size(); i++)
            sum += cs[i].numberOfVariants();
        return sum;
    }
};

class Store
{
private:
    ptree c_old;
    int lineCount;
    // column width for output; how often (every x lines) to print table header
    enum {columnWidth=14, tableHeaderEveryXLines=20};

    typedef std::pair<std::string, double> RPair;
    typedef std::vector<RPair> RVector;

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
        std::vector<std::string> names = getVParams(c);
        for (int i=0; i<names.size(); i++)
        {
            ptree::assoc_iterator j = c.find(names[i]);
            if (j != c.not_found())
                c.erase(c.to_iterator(j));
        }
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
        write_json(ss, c);
        std::string s = ss.str();
        prependLines(s, "# ");
        std::cout << s;
        std::cout << "# " << std::endl;
        
        lineCount = 0;
    }

    void printTableHeader (const ptree& c, const rtree& r) const
    {
        std::vector<std::string> ps = getVParams(c);
        RVector os = getFlattenedTree(r);
        std::cout << "# ";
        for (int i=0; i<ps.size(); i++)
            std::cout << std::setw(columnWidth) << ps[i];
        for (int i=0; i<os.size(); i++)
            std::cout << std::setw(columnWidth) << os[i].first;
        std::cout << std::endl;
    }

    void printTableRow (const ptree& c, const rtree& r)
    {
        if (lineCount % tableHeaderEveryXLines == 0)
            printTableHeader(c,r);

        std::vector<std::string> ps = getVParams(c);
        RVector os = getFlattenedTree(r);
        std::cout << "  ";
        for (int i=0; i<ps.size(); i++)
            std::cout << std::setw(columnWidth) << c.get<double>(ps[i], 0);
        for (int i=0; i<os.size(); i++)
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

    RVector getFlattenedTree (const rtree& r) const
    {
        RVector v;
        getFlattenedTreeRecursive(r, v, "");
        return v;
    }
    
    void getFlattenedTreeRecursive (const rtree& r, RVector& v, const std::string& path) const
    {
        if (r.size() == 0)
        {
            v.push_back(RPair(path, r.data()));
            return;
        }

        for (rtree::const_iterator i=r.begin(); i!=r.end(); i++)
            if (path.length() == 0)
                getFlattenedTreeRecursive(i->second, v, i->first);
            else
                getFlattenedTreeRecursive(i->second, v, path + "." + i->first);
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

    // std::stringstream prependLines (const std::stringstream& is, const std::string& pr) const
    // {
    //     std::stringstream os;
    //     std::string s;
    //     while (!is.eof())
    //     {
    //         os << pr;
    //         std::getline(is, s);
    //         os << s;
    //     }
    //     return os;
    //     // std::string s;
    //     // ss.seekg(0);
    //     // ss.seekp(0);
    //     // while (!ss.eof())
    //     // {
    //     //     ss << pr;
    //     //     std::getline(ss, s);
    //     //     ss.seekp(ss.tellg());
    //     // }
    //     // ss.seekg(0);
    // }
};

#endif // __CQCA_HPP__
