/**
 * @file 
 * Configurable QCA classes and assorted functionality.
 *
 * Boost's property tree (boost::property_tree::ptree) is used to store
 * configuration data. In the code, "ptree" and "configuration" are often used
 * interchangeably. A result tree (rtree), a variant of ptree, is used to store
 * computed results.
 *
 * The file contains: 
 * * PropertyTree (useful helper funtions for ptree)
 * * the configurable QCA classes (CQca, CLayout, etc)
 * * a system to generate multiple QCA system configurations out of one
 *   single configuration file (Configurator, VConfiguration, VParam)
 * * a Store, for nicely formatted output
 * * a stochastic optimization method to find maxima (StochasticFindMax, Random)
 * * Runner, to run the simulation
 *
 * @author Burkhard Ritter <burkhard@ualberta.ca>
 */
#ifndef __CQCA_HPP__
#define __CQCA_HPP__

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <ctime>
#include "qca.hpp"
#include "version.hpp"

class ConfigurationException: public std::runtime_error
{
public:
    ConfigurationException(const std::string& msg)
        : std::runtime_error(msg) {}
};

/**
 * Functionality related to boost::property_tree.
 *
 * Array convention
 * ----------------
 *
 * I use a different convention for representing arrays in a ptree than
 * boost::property_tree. In boost::property_tree an array is a ptree node where
 * all children have an empty key "". I use the convention where an array is a
 * ptree node where the keys of all children are enumerated, "0", "1", "2", and
 * so on; and the enumeration must match the order of the children. This
 * convention makes array element access much easier. Hence on reading and
 * writing ptrees, I convert between the two conventions.
 * @see treeFromJson
 * @see treeToJson
 */
namespace PropertyTree
{

using boost::property_tree::ptree;

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

/**
 * Get the children of the given ptree as an array.
 *
 * Not recursive, only direct children.
 */
template<typename T>
std::vector<T> getArray (const ptree& c)
{
    std::vector<T> v;
    for (ptree::const_iterator i=c.begin(); i!=c.end(); i++)
        v.push_back(i->second.get_value<T>());
    return v;
}

/**
 * Construct a ptree (or rtree) representing the given array.
 *
 * Uses my ptree array convention.
 * @see PropertyTree
 */
template <class Tree, typename T>
Tree constructArray (const std::vector<T>& v)
{
    typedef typename Tree::value_type pair;
    Tree c;
    for (size_t i=0; i<v.size(); i++)
    {
        Tree p;
        p.put_value(v[i]);
        c.push_back(pair(toString(i), p));
    }
    return c;
}

/**
 * Check whether the given ptree (or rtree) represents an array.
 *
 * Uses my ptree array convention.
 * @see PropertyTree
 */
template<class Tree>
bool isArray (const Tree& c)
{
    if (c.data() != typename Tree::data_type())
        return false;
    if (c.size() == 0)
        return false;
    bool flag = true;
    int j=0;
    for (typename Tree::const_iterator i=c.begin(); i!=c.end(); i++)
    {
        if (i->first != toString(j)) flag = false;
        j++;
    }
    return flag;
}

/**
 * Does the given tree node has a direct child with the given key.
 */
template<class Tree>
bool hasKey (const Tree& t, const typename Tree::key_type& k)
{
    typename Tree::const_assoc_iterator i = t.find(k);
    return i != t.not_found();
}

namespace detail
{

/**
 * Function object that converts a ptree node from the
 * boost::property_tree array convention to my array convention.
 * @see PropertyTree
 */
struct ConvertArray
{
    void operator() (ptree& c)
    {
        // is it an array
        if (c.data() != "")
            return;
        if (c.size() == 0)
            return;
        bool flag = true;
        for (ptree::const_iterator i=c.begin(); i!=c.end(); i++)
            if (i->first != "") flag = false;

        // change indices from "" to a proper index
        if (flag)
        {
            std::vector<ptree> cs;
            for (ptree::const_iterator i=c.begin(); i!=c.end(); i++)
                cs.push_back(i->second);
            c.erase(c.begin(), c.end());
            for (size_t i=0; i<cs.size(); i++)
                c.put_child(toString(i), cs[i]);
        }
    }
};

/**
 * Function object that converts a ptree node from my array convention
 * to the boost::property_tree array convention.
 * @see PropertyTree
 */
struct ConvertArrayBack
{
    void operator() (ptree& c)
    {
        // is it an array
        if (c.data() != "")
            return;
        if (c.size() == 0)
            return;
        bool flag = true;
        int j = 0;
        for (ptree::const_iterator i=c.begin(); i!=c.end(); i++)
        {
            if (i->first != toString(j)) flag = false;
            j++;
        }

        // change indices from numbers to ""
        if (flag)
        {
            std::vector<ptree> cs;
            for (ptree::const_iterator i=c.begin(); i!=c.end(); i++)
                cs.push_back(i->second);
            c.erase(c.begin(), c.end());
            for (size_t i=0; i<cs.size(); i++)
                c.push_back(ptree::value_type("", cs[i]));
        }
    }
};

/**
 * Walk the ptree recursively and call the function object for each
 * node.
 */
template<class Function>
void walkTree (ptree& c, Function f)
{
    for (ptree::iterator i=c.begin(); i!=c.end(); i++)
    {
        walkTree(i->second, f);
    }
    f(c);
}

/**
 * Converts "simplified json" to proper json.
 *
 * Simplified json makes quoting ("") keys and values optional by
 * treating all keys and values as json strings.
 */
std::string jsonify (std::istream& is)
{
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

std::string jsonify (const std::string& s)
{
    std::stringstream is(s);
    return jsonify(is);
}

bool isfile (std::string path)
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

} // namespace detail

/**
 * Read a property tree from "simplified json".
 *
 * Simplified json makes quoting ("") keys and values optional by
 * treating all keys and values as json strings.
 * @see jsonify
 *
 * Converts the ptree from the boost::property_tree array convention to my
 * convention.
 * @see PropertyTree
 *
 * @param json Either a simplified json string or a filename of a file
 * containing simplified json.
 * @return The property tree.
 */
ptree treeFromJson (const std::string& json)
{
    ptree c;
    if (detail::isfile(json))
    {
        std::ifstream is(json.c_str());
        std::stringstream ss(detail::jsonify(is));
        is.close();
        read_json(ss, c);
    }
    else
    {
        std::stringstream ss(detail::jsonify(json));
        read_json(ss, c);
    }
    detail::walkTree(c, detail::ConvertArray());
    return c;
}

/**
 * Convert the property tree to a json string.
 *
 * Converts the ptree from my array convention to the boost::property_tree
 * convention.
 * @see PropertyTree
 */
std::string treeToJson (ptree c)
{
    std::stringstream ss;
    detail::walkTree(c, detail::ConvertArrayBack());
    write_json(ss, c);
    return ss.str();
}

} // namespace PropertyTree

/**
 * Result tree.
 *
 * Used to store results of computations. Similar to
 * boost::property_tree::ptree, but all values are of type double instead of
 * type string.
 */
typedef boost::property_tree::basic_ptree<std::string, double> rtree;

namespace PT = PropertyTree;
using boost::property_tree::ptree;

/**
 * Configurable layout.
 */
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
        if (PT::hasKey(c, "a") && PT::hasKey(c, "V1"))
            throw ConfigurationException("Either specify 'a' or 'V1', but "
                                         "not both at the same time");
        if (PT::hasKey(c, "V1"))
            return 1.0 / c.get<double>("V1");
        return c.get<double>("a", 1.0);
    }
    
    double getB (const ptree& c) const
    {
        if (PT::hasKey(c, "b") && PT::hasKey(c, "boa"))
            throw ConfigurationException("Either specify 'b' or 'boa', but "
                                         "not both at the same time");
        if (PT::hasKey(c, "boa"))
        {
            double a_ = getA(c);
            return c.get<double>("boa") * a_;
        }
        return c.get<double>("b", 3.0);
    }
    
    std::vector<double> getBs (const ptree& c) const
    {
        if (PT::hasKey(c, "bs") && PT::hasKey(c, "boas"))
            throw ConfigurationException("Either specify 'bs' or 'boas', but "
                                         "not both at the same time");
        if (PT::hasKey(c, "boas"))
        {
            double a_ = getA(c);
            std::vector<double> bs_ = PT::getArray<double>(c.get_child("boas"));
            for (size_t i=0; i<bs_.size(); i++)
                bs_[i] *= a_;
            return bs_;
        }
        return PT::getArray<double>(c.get_child("bs", ptree()));
    }

    void setAOrV1 (ptree& c) const
    {
        if (PT::hasKey(c, "V1"))
            c.put("V1", 1.0 / a);
        else
            c.put("a", a);
    }
    
    void setBOrBoa (ptree& c) const
    {
        if (PT::hasKey(c, "boa"))
            c.put("boa", b / a);
        else
            c.put("b", b);
    }

    void setBsOrBoas (ptree& c) const
    {
        if (PT::hasKey(c, "boas"))
        {
            std::vector<double> boas = bs;
            for (size_t i=0; i<boas.size(); i++)
                boas[i] /= a;
            c.put_child("boas", PT::constructArray<ptree>(boas));
        }
        else
            c.put_child("bs", PT::constructArray<ptree>(bs));
    }
};

/**
 * Configurable QCA base class.
 */
class CQcaGenericBase
{
public:
    virtual ~CQcaGenericBase () {}
    virtual void setConfig (const ptree& c) = 0;
    virtual ptree getConfig () const = 0;
    virtual rtree measure () = 0;
};

/**
 * Generic configurable QCA class.
 *
 * The template argument is used to specify which model to use (that is, the specific
 * QCA class, e.g. QcaBond, QcaFixedCharge, or QcaGrandCanonical).
 */
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
            else if (i->first == "DOS")
                r.put_child("DOS", measureDOS(i->second));
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
        else if (PT::isArray(c))
            cells = PT::getArray<int>(c);
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
                        PT::toString(cells[i]) + ", but there are only " + 
                        PT::toString(s.N_p()) + " cells in the system.");
        }

        return cells;
    }

    rtree measureP (const ptree& c) const
    {
        std::vector<int> cells = getCells(c);
        rtree r;
        for (size_t i=0; i<cells.size(); i++)
            r.put(PT::toString(cells[i]), s.measurePolarization(cells[i]));
        return r;
    }

    rtree measureP2 (const ptree& c) const
    {
        std::vector<int> cells = getCells(c);
        rtree r;
        for (size_t i=0; i<cells.size(); i++)
            r.put(PT::toString(cells[i]), s.measurePolarization2(cells[i]));
        return r;
    }

    rtree measureN (const ptree& c) const
    {
        std::vector<int> cells = getCells(c);
        rtree r;
        for (size_t i=0; i<cells.size(); i++)
        {
            std::vector<double> ns = s.measureParticleNumber(cells[i]);
            r.put(PT::toString(cells[i]) + ".total", ns[4]);
            r.put(PT::toString(cells[i]) + ".0", ns[0]);
            r.put(PT::toString(cells[i]) + ".1", ns[1]);
            r.put(PT::toString(cells[i]) + ".2", ns[2]);
            r.put(PT::toString(cells[i]) + ".3", ns[3]);
        }
        return r;
    }

    rtree measureE (const ptree& c) 
    {
        rtree r;
        for (int i=0; i<s.energies().size(); i++)
        {
            rtree v;
            v.put("abs", s.energies()(i));
            v.put("rel", s.energies()(i) - s.Emin());
            r.push_back(rtree::value_type(PT::toString(i), v));
        }
        return r;
    }

    rtree measureDOS (const ptree& c)
    {
        rtree r;
        double deltaE = c.get<double>("deltaE", 0.1);
        double Emin = c.get<double>("Emin", s.energies().minCoeff());
        double Emax = c.get<double>("Emax", s.energies().maxCoeff());
        int Nmin = Emin / deltaE;
        int Nmax = Emax / deltaE;
        if (Nmin * deltaE > Emin) Nmin--;
        if (Nmax * deltaE < Emax) Nmax++;
        Emin = Nmin * deltaE;
        Emax = Nmax * deltaE;

        std::vector<int> h(Nmax-Nmin+1, 0);
        for (int i=0; i<s.energies().size(); i++)
        {
            if (s.energies()(i) >= Emin && s.energies()(i) <= Emax)
                h[(s.energies()(i) - Emin)/deltaE]++;
            else
                std::cerr 
                    << "Warning: Energy value outside of DOS energy range." 
                    << std::endl;
        }
        for (size_t i=0; i<h.size(); i++)
        {
            rtree n;
            n.put("E", i*deltaE + Emin);
            n.put("DOS", static_cast<double>(h[i])/s.energies().size());
            r.push_back(rtree::value_type(PT::toString(i), n));
        }
        return r;
    }

    double getBeta (const ptree& c) const
    {
        if (PT::hasKey(c, "beta") && PT::hasKey(c, "T"))
            throw ConfigurationException("Either specify 'beta' or 'T', but "
                                         "not both at the same time");
        if (PT::hasKey(c, "T"))
            return 1.0 / c.get<double>("T");
        return c.get<double>("beta", 1);
    }

    void setBetaOrT (ptree& c, double beta) const
    {
        if (PT::hasKey(c, "T"))
            c.put("T", 1.0 / beta);
        else
            c.put("beta", beta);
    }
};

/**
 * Configurable QCA class.
 *
 * The class takes a configuration describing any valid qca system and what to
 * calculate/measure via setConfig(). It instantiates the correct model (bond,
 * fixed, grand canonical) and configures the qca system correctly. Then the
 * observables of intererst (as specified in the configuration) can be obtained
 * as a result tree by calling measure().
 * @see rtree
 */
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

/**
 * Varying parameter.
 *
 * @see Configurator
 *
 * A VParam is either a list (json array) or a list generated from a
 * "generator". The generator syntax is start-value:end-value:interval. 
 *
 * increment() increments to the next element in the list. value() gets the
 * current element's value. name() is the name of the VParam. size() returns the
 * number of elements in the list.
 */
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
    size_t index;

public:
    VParam (ptree&  c, const std::string& name__)
    : name_(name__)
    {
        try {
            c.get_child(name_);
        }
        catch (boost::property_tree::ptree_bad_path e)
        {
            throw ConfigurationException("Changing parameter '" + name_ + 
                                         "' does not exist.");
        }

        ptree& p = c.get_child(name_);
        if (PT::isArray(p))
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
            n.put_value(PT::toString(values[index]));
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
            T item = PT::fromString<T>(v.substr(pos1, pos2-pos1));
            list.push_back(item);
            pos1 = pos2+1;
        }
        return list;
    }
};

/**
 * Configuration with Variants.
 *
 * @see Configurator
 *
 * Each VConfiguration can have multiple "Variants". A VConfiguration has
 * varying parameters (VParam), as specified in the "changing" json array. A
 * VParam is effectively a parameter with a list of different values. A
 * "Variant" of a VConfiguration corresponds to one specific set of values of
 * all possible sets of values of VParams.
 */
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
        if (!PT::isArray(p))
        {
            if (p.get_value<std::string>() != "")
                throw ConfigurationException("'changing' must be an array of changing parameters");
            hasVariants = false;
            return;
        }

        hasVariants = true;
        changed.clear();
        std::vector<std::string> names = PT::getArray<std::string>(p);
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
        nc.put_child("changed", PT::constructArray<ptree>(changed));
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

/**
 * Generates individual qca system configurations from one single global
 * configuration file or string.
 *
 * The single global configuration is a list of qca system configurations,
 * possibly using parameter lists or generators. Individual qca system
 * configurations (which can be fed to the CQca class) can be retrieved with
 * getNext().
 *
 * More specifically, the single global configuration is either a json list (an
 * array) of configurations or just one configuration by itself. Each of these
 * configurations (represented by VConfiguration) can have "Variants". Each
 * Variant then is one final concrete qca system configuration. A VConfiguration
 * has one or more varying parameters (represented by VParam). In JSON these
 * parameters must be listed in the "changing" array. Each VParam is either a
 * list of values (a JSON array) or a "generator" that generates this list. The
 * generator syntax is start-value:end-value:interval. See the test cases in
 * cqcaTest.cpp for examples.
 */
class Configurator 
{
private:
    std::vector<VConfiguration> cs;
    size_t i_c;

public:
    Configurator (const std::string& json)
    : i_c(0)
    {
        ptree c = PT::treeFromJson(json);
        if (PT::isArray(c))
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
};

/**
 * Write results of measurements to the standard output.
 *
 * "Stores" measurements to stdout. Results of measurements are held in a result
 * tree (rtree). The configuraion that was used to calculate those results is
 * held in a property tree (ptree). Prints nicely formatted output and also
 * understands varying parameters (VParam) and Variants of configurations
 * (VConfiguration). 
 */
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
        // remove "internal" fields
        c.put_child("changed", ptree());
        c.put_child("_additionalcolumns", ptree());
        return c;
    }

    void printHeader (const ptree& c)
    {
        // only false if this is the first header / measurement
        if (lineCount != 0)
            std::cout << std::endl << std::endl;

        std::cout << "# program version: " << GIT_PROGRAM_VERSION << std::endl 
                  << "# date: " << getDate() << std::endl
                  << "# " << std::endl;
        
        std::stringstream ss(PT::treeToJson(prepareForPrinting(c)));
        std::string line;
        while (getline(ss, line))
        {
            std::cout << "# " << line << std::endl;
        }
        std::cout << "# " << std::endl;
        
        lineCount = 0;
    }

    void printTableHeader (const ptree& c, const rtree& r) const
    {
        std::cout << "# ";

        if (PT::hasKey(c, "_additionalcolumns"))
        {
            const ptree ac = c.get_child("_additionalcolumns");
            for (ptree::const_iterator i=ac.begin(); i!=ac.end(); i++)
                std::cout << std::setw(columnWidth) << i->first;
        }
        
        std::vector<std::string> ps = getVParams(c);
        for (size_t i=0; i<ps.size(); i++)
        {
            const ptree& n = c.get_child(ps[i]);
            if (PT::isArray(n))
                for (size_t j=0; j<n.size(); j++)
                    std::cout << std::setw(columnWidth) << (ps[i] + "." + PT::toString(j));
            else 
                std::cout << std::setw(columnWidth) << ps[i];
        }
        
        TreeTypes<rtree>::vector os = getFlattenedTree(r);
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

        std::cout << "  ";
        
        if (PT::hasKey(c, "_additionalcolumns"))
        {
            const ptree ac = c.get_child("_additionalcolumns");
            for (ptree::const_iterator i=ac.begin(); i!=ac.end(); i++)
                std::cout << std::setw(columnWidth) << i->second.get_value<double>();
        }
        
        std::vector<std::string> ps = getVParams(c);
        for (size_t i=0; i<ps.size(); i++)
        {
            const ptree& n = c.get_child(ps[i]);
            if (PT::isArray(n))
                for (ptree::const_iterator j=n.begin(); j!=n.end(); j++)
                    std::cout << std::setw(columnWidth) << j->second.get_value<double>();
            else
                std::cout << std::setw(columnWidth) << c.get<double>(ps[i], 0);
        }
        
        TreeTypes<rtree>::vector os = getFlattenedTree(r);
        for (size_t i=0; i<os.size(); i++)
            std::cout << std::setw(columnWidth) << os[i].second;
        
        std::cout << std::endl;
        lineCount++;
    }

    std::vector<std::string> getVParams (const ptree& c) const
    {
        std::vector<std::string> vparams;
        ptree p1 = c.get_child("changing", ptree());
        if (PT::isArray(p1))
        {
            const std::vector<std::string> a1 = PT::getArray<std::string>(p1);
            vparams.insert(vparams.end(), a1.begin(), a1.end());
        }
        ptree p2 = c.get_child("findmax.optimize", ptree());
        if (PT::isArray(p2))
        {
            const std::vector<std::string> a2 = PT::getArray<std::string>(p2);
            vparams.insert(vparams.end(), a2.begin(), a2.end());
        }
        return vparams;
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

        bool isarray = PT::isArray(t);
        int l=0;
        for (typename Tree::const_iterator i=t.begin(); i!=t.end(); i++)
        {
            std::string npath = path;
            if (npath.length() != 0)
                npath += ".";
            if (isarray)
                npath += PT::toString(l);
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
    
    ptree prepareForPrinting (ptree c)
    {
        ptree::assoc_iterator i;
        i = c.find("_additionalcolumns");
        if (i != c.not_found())
            c.erase(c.to_iterator(i));
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

class Random
{
private:
    typedef boost::mt19937 baseGeneratorT;
    typedef boost::uniform_real<> distributionT;
    typedef boost::variate_generator<baseGeneratorT&, distributionT> generatorT;

    baseGeneratorT baseGenerator;
    distributionT distribution;
    generatorT generator;

public:
    Random (unsigned int seed=42u)
    : baseGenerator(seed), distribution(0,1), 
      generator(baseGenerator, distribution)
    {
    }
    
    double operator() ()
    {
        return generator();
    }
};

/**
 * Stochastic optimization method to find maxima.
 *
 * This is a stochastic annealing optimization method to find the maximum of a
 * chosen observable, varying a chosen set of parameters. E.g. Optimize the
 * inter-cell spacings in a non-uniform wire for maximal output polarization.
 *
 * The method is from the following paper:
 * * PhysRevB.76.104432
 *    * title: Variational ground states of two-dimensional antiferromagnets in the valence bond basis
 *    * author: Lou, Jie and Sandvik, Anders W.
 *    * journal: Phys. Rev. B
 *    * volume: 76
 *    * issue: 10
 *    * pages: 104432
 *    * numpages: 8
 *    * year: 2007
 *    * month: Sep
 *    * doi: 10.1103/PhysRevB.76.104432
 *    * url: http://link.aps.org/doi/10.1103/PhysRevB.76.104432
 *    * publisher: American Physical Society
 *
 * Usage
 * -----
 *
 * In the json configuration the following fields must be defined:
 * * findmax.of -- find maximum of this observable, e.g. P.2
 * * findmax.optimize -- list of parameters that should be varied / optimized,
 *                       e.g. layout.boas
 *
 * The following fields are optional:
 * * findmax.iterations -- number of iteration steps, default is 1000
 */
class StochasticFindMax
{
private:
    Random R;
    
    // constant after initialization
    std::string findMaxOf;
    std::vector<std::string> ans; // argument names
    int maxN;
    double alpha;

    // change with each iteration step
    double beta;
    std::vector<double> avs; // current argument values
    std::vector<double> oldAvs; // old argument values
    double fv; // current function value
    double oldFv; // old function value
    
    ptree c;
    CQca& s;
    Store& o;
    rtree r;

public:
    StochasticFindMax (ptree c_, CQca& s_, Store& o_)
    : alpha(0.75), beta(1), c(c_), s(s_), o(o_)
    {
        readFindMaxConfig();

        // initialize old argument values, function value
        oldAvs.resize(avs.size());
        for (size_t i=0; i<avs.size(); i++)
            oldAvs[i] = avs[i] + 0.01 * avs[i];
        calculateFunctionValue();
        oldFv = fv;
    }

    /**
     * Find the maximum.
     *
     * Implements the annealing loop and uses the following annealing scheme:
     *
     * @f[ \beta = \frac{1}{i^\alpha} @f]
     * 
     * with i the iteration number 
     * and @f$ 0.5 < \alpha < 1 @f$, 
     * @f$ \alpha \sim 0.75 @f$ is recommended as a starting point
     */
    void findmax ()
    {
        for (size_t i=1; i<=maxN; i++)
        {
            beta = std::pow(i, -alpha);

            // additional columns in output
            c.put("_additionalcolumns.Iteration", i);
            c.put("_additionalcolumns.beta", beta);
            
            update();
            o.store(s.getConfig(), r);
        }
    }

private:
    /**
     * Calculates the function value (findmax.of) for the current argument
     * values (findmax.optimize). 
     *
     * Updates the variables fv and r.
     */
    void calculateFunctionValue ()
    {
        assert (ans.size() == avs.size());
        for (size_t i=0; i<ans.size(); i++)
            c.put(ans[i], avs[i]);
        s.setConfig(c);
        r = s.measure();
        fv = r.get<double>(findMaxOf);
    }

    /**
     * Update the argument values and the function value.
     *
     * Implements the stochastic scheme as detailed in the paper. 
     * @see StochasticFindMax
     *
     * More specifically, the update scheme is, for each argument
     * @f[ \ln av^{i+1} = \ln av^i + R \beta \mathrm{sign} 
     *     \left( \frac{\partial fv}{\partial av} \right) @f]
     * where R is a random number [0,1] and @f$ \beta @f$ is the annealing
     * temperature.
     *
     */
    void update ()
    {
        oldFv = fv;

        // determine sign of slope for all arguments
        std::vector<int> ss(avs.size());
        for (size_t i=0; i<ss.size(); i++)
        {
            const double ov = avs[i];
            avs[i] += 0.1 * std::abs(avs[i] - oldAvs[i]);
            calculateFunctionValue();
            avs[i] = ov;
            if (fv >= oldFv)
                ss[i] = 1;
            else
                ss[i] = -1;
        }

        // update argument values
        oldAvs = avs;
        for (size_t i=0; i<ss.size(); i++)
            avs[i] = avs[i] * std::exp(R() * beta * ss[i]);

        // update function value
        calculateFunctionValue();
    }

    void readFindMaxConfig ()
    {
        findMaxOf = c.get<std::string>("findmax.of");
        maxN = c.get("findmax.iterations", 1000);

        // which arguments to optimize ("findmax.optimize")
        // read their names to ans and their initial values to avs
        ans.clear();
        avs.clear();
        const ptree& p = c.get_child("findmax.optimize");
        for (ptree::const_iterator i=p.begin(); i!=p.end(); i++)
        {
            const std::string name = i->second.get_value<std::string>();
            // "unfold" arrays, i.e. [a] -> a.0, a.1, a.2 and so on
            if (PT::isArray(c.get_child(name)))
            {
                std::vector<double> vs = PT::getArray<double>(c.get_child(name));
                for (size_t j=0; j<vs.size(); j++)
                {
                    ans.push_back(name + "." + PT::toString(j));
                    avs.push_back(vs[j]);
                }
            }
            else
            {
                ans.push_back(name);
                avs.push_back(c.get<double>(name));
            }
            ptree op = c.get_child(name);
            c.put_child("original." + name, op);
        }
    }
};

/**
 * Run the simulation.
 *
 * Uses Configurator to retrieve all qca system configurations (including
 * Variants), CQca to do the actual calculations and Store to print the results
 * to the standard output.
 */
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
            const ptree cc = c.getNext();
            
            if (PT::hasKey(cc, "findmax"))
            {
                StochasticFindMax f(cc, s, o);
                f.findmax();
            }
            else
            {
                s.setConfig(cc);
                o.store(s.getConfig(), s.measure());
            }
        }
    }
};

#endif // __CQCA_HPP__
