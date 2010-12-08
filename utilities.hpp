/*
 * utilities.hpp
 *
 * Originally written for the ctqmc project.
 *
 *  Created on: Jul 31, 2009
 *  Author: "Burkhard Ritter"
 */

#ifndef __UTILITIES_HPP__
#define __UTILITIES_HPP__

#include <map>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <cctype>
namespace System {
    #include <sys/types.h>
    #include <sys/stat.h>
}

/**
 * The interface follows the python os.path module.
 */
namespace FileSystem
{
    namespace Path
    {
        inline bool exists(std::string path)
        {
            struct System::stat s;
            if (System::stat(path.c_str(), &s) == 0) {
                return true;
            }
            return false;
        }

        inline bool isdir(std::string path)
        {
            struct System::stat s;
            if (System::stat(path.c_str(), &s) != 0) {
                return false;
            }
            if (S_ISDIR(s.st_mode)) {
                return true;
            }
            return false;
        }

        inline bool isfile(std::string path)
        {
            struct System::stat s;
            if (System::stat(path.c_str(), &s) != 0) {
                return false;
            }
            if (S_ISREG(s.st_mode)) {
                return true;
            }
            return false;
        }

    }

    class FileSystemException: public std::runtime_error
    {
    public:
        FileSystemException(const std::string& msg)
            : std::runtime_error(msg) {}
    };

    inline void mkdir(std::string path)
    {
        if (System::mkdir(path.c_str(), 0755) != 0) {
            throw FileSystemException(
                std::string("Could not create directory '") + path + "'."
            );
        }
    }
}

namespace Helpers
{

inline std::string indent(std::string toIndent, std::string indentation)
{
    for (size_t pos=0; pos<toIndent.length(); ) {
        toIndent.insert(pos, indentation);
        pos += indentation.length();
        pos = toIndent.find('\n', pos);
        if (pos == std::string::npos) break;
        pos++;
    }
    return toIndent;
}

inline std::string trim(std::string s)
{
    std::string whitespaces (" \t\f\v\n\r");
    size_t pos1 = s.find_first_not_of(whitespaces);
    size_t pos2 = s.find_last_not_of(whitespaces);
    if (pos1 == std::string::npos || pos2 == std::string::npos) {
        //s consists only of whitespaces
        s.clear();
    }
    else {
        s = s.substr(pos1, pos2-pos1+1);
    }
    return s;
}

inline std::string toLower (std::string s)
{
    std::transform(s.begin(), s.end(), s.begin(), tolower);
    return s;
}

}

using namespace Helpers;

class ConversionException : public std::exception {};

/*
 * TODO: more detailed error messages
 */
class DescriptionItem
{
public:
    DescriptionItem(): value(), set(false) {}
    DescriptionItem(std::string value_, bool set_ = true): value(value_), set(set_) {}
    DescriptionItem(int value_, bool set_ = true): set(set_) {setValue(value_);}
    DescriptionItem(size_t value_, bool set_ = true): set(set_) {setValue(value_);}
    DescriptionItem(double value_, bool set_ = true): set(set_) {setValue(value_);}

    virtual ~DescriptionItem() {}

    template<typename T>
    void operator= (const T& v)
    {
        setValue(v);
        set = true;
    }

    operator double() const
    {
        double returnValue;
        std::stringstream s(value);
        s >> returnValue;
        if (!s.eof() || value.empty()) {
            throw ConversionException();
        }
        return returnValue;
    }

    operator int() const
    {
        /*
        //both, sscanf as well as as atoi don't work reliably
        //they don't give errors when you try to convert, say, 2.4 to int
        //in this case they just return 2
        int returnValue;
        if (std::sscanf(value.c_str(), "%d", &returnValue) == 1) {
            return returnValue;
        }

        //this is using boost::lexical_cast and works very well, but we get the
        //boost dependency
        try {
            return boost::lexical_cast<int>(value);
        }
        catch(boost::bad_lexical_cast& e) {
            throw ConversionException();
        }
        */

        //this works and avoids the boost dependency
        int returnValue;
        std::stringstream s(value);
        s >> returnValue;
        //check if the stream is empty - it should be if it only contained an
        //integer
        if (!s.eof() || value.empty()) {
            throw ConversionException();
        }
        return returnValue;
    }

    operator size_t() const
    {
        size_t returnValue;
        std::stringstream s(value);
        s >> returnValue;
        if (!s.eof() || value.empty()) {
            throw ConversionException();
        }
        return returnValue;
    }

    /*
    operator unsigned int()
    {
        uint32_t returnValue;
        std::stringstream s(value);
        s >> returnValue;
        if (!s.eof()) {
            throw ConversionException();
        }
        return returnValue;
    }
    */

    operator bool() const
    {
        std::string lvalue = toLower(value);
        if (lvalue == "no" || lvalue == "false")
            return false;
        else if (lvalue == "yes" || lvalue == "true")
            return true;

        throw ConversionException();
    }

    operator std::string() const
    {
        return value;
    }
    
    /*
     * TODO: obviously, std::vector won't scale up indefinitely
     */
    /*
     * TODO: would be nice to implement an alternative syntax
     * 1,2,3,4,5 for a simple list
     * 1:2:0.2 for a generated list
     * hence we would not need escaping on the shell
     */
    template<typename T>
    operator std::vector<T> () const
    {
        if (value.size() == 0)
            return std::vector<T>(); //throw ConversionException();

        if (value[0] == '[')
        {
            if (value.find_first_of("]") != value.size()-1)
                throw ConversionException();
            return getList<T>(value.substr(1, value.size()-2));
        }
        if (value[0] == '(')
        {
            if (value.find_first_of(")") != value.size()-1)
                throw ConversionException();
            std::vector<T> tuple = getList<T>(value.substr(1, value.size()-2));
            if (tuple.size() != 3) throw ConversionException();
            T begin = tuple[0];
            T end = tuple[1];
            T inc = tuple[2];
            if (inc==0 || (begin<=end && inc<0) || (end<=begin && inc>0)) 
                throw ConversionException();
            std::vector<T> generated;
            for (T i=begin; (begin<=end && i<=end) || (end<=begin && i>=end); i+=inc)
                generated.push_back(i);
            return generated;
        }
        else
            return std::vector<T>(1, get<T>());
    }

    bool operator== (const std::string& s) const
    {
        return value == s;
    }

    bool operator!= (const std::string& s) const
    {
        return value != s;
    }

    /*
     * This is only an alias for direct assignments.
     *
     * T blub = instanceOfIniValue;
     * and
     * T blub = instanceOfIniValue.get();
     * should be equivalent.
     */
    template<typename T>
    T get() const
    {
        return *this;
    }

    template<typename T>
    T get(const T& defaultValue) const
    {
        if (!set) return defaultValue;
        return *this;
    }

    /**
     * Set a default value, that is, set a value, if this DescriptionItem doesn't
     * have a value yet.
     *
     * @param defaultValue
     */
    template<typename T>
    void setDefault(const T& defaultValue)
    {
        if (!set) setValue(defaultValue);
    }

    bool isSet () const
    {
        return set;
    }

    const std::string& getValue() const {return value;}

private:
    template<typename T>
    void setValue(const T& v)
    {
        std::stringstream s;
        s << v;
        value = trim(s.str());
    }

    //TODO: could we also use template specialisation?
    //template<>
    //void setValue<bool>(const bool& v)
    void setValue(const bool& v)
    {
        if (v) value = "true";
        else value = "false";
    }

    template<typename T>
    std::vector<T> getList (const std::string& v) const
    {
        if (v.size() == 0) return std::vector<T>();
        std::vector<T> list;
        DescriptionItem item;
        size_t pos1=0;
        size_t pos2=0;
        while (pos1 <= v.size())
        {
            pos2 = v.find_first_of(',', pos1);
            if (pos2 == std::string::npos) pos2 = v.size();
            if (pos2 == pos1) throw ConversionException();
            item = v.substr(pos1, pos2-pos1);
            list.push_back(item);
            pos1 = pos2+1;
        }
        return list;
    }

    std::string value;
    bool set;
};

inline std::ostream& operator<< (std::ostream& o, const DescriptionItem& v)
{
    o << v.getValue();
    return o;
}

class Description
{
public:
    typedef std::map<std::string, Description> SectionsType;
    typedef std::map<std::string, DescriptionItem> ItemsType;

    Description(std::string name)
        : name(name) {}

    //needed to store Description objects in a map
    Description() {}

    virtual ~Description() {}

    bool hasSection(std::string sectionName)
    {
        return sections.count(sectionName);
    }

    Description getSection(std::string sectionName)
    {
        if (!hasSection(sectionName)) {
            sections[sectionName] = Description(sectionName);
        }
        return sections[sectionName];
    }

    bool hasItem(std::string itemName)
    {
        return items.count(itemName);
    }

    DescriptionItem getItem(std::string itemName)
    {
        return items[itemName];
    }

    DescriptionItem& operator[] (std::string itemName)
    {
        return items[itemName];
    }

    std::string asString() const
    {
        std::stringstream s;
        s << "[" << name << "]" << std::endl;
        typedef ItemsType::const_iterator IIT;
        for (IIT i=items.begin(); i!=items.end(); i++) {
            s << i->first << "\t= " << i->second << std::endl;
        }
        if (items.size() > 0) {
            s << std::endl;
        }

        typedef SectionsType::const_iterator SIT;
        for (SIT i=sections.begin(); i!=sections.end(); i++) {
            s << indent(i->second.asString(), "  ") << std::endl;
        }

        return s.str();
    }

    std::string name;
    SectionsType sections;
    ItemsType items;
};

class IniException: public std::runtime_error
{
public:
    IniException(const std::string message)
        : std::runtime_error(message) {}
};

class IniFile: public Description
{
public:
    //convenience constructor
    IniFile()
        : Description("IniFile") {}

    IniFile(std::string fileName)
        : Description("IniFile"), fileName(fileName)
    {
        parseIniFile();
    }

    virtual ~IniFile() {}

private:
    void parseIniFile()
    {
        std::ifstream file(fileName.c_str());
        if(!file) {
            throw IniException("Configuration file '" + fileName + "' does not exist.");
        }

        std::string line, sectionName, parentSectionName, key, value;
        while (std::getline(file, line)) {
            /*
             * skip comments
             */
            if (line[0] == '#') {
                continue;
            }
            /*
             * a section
             *
             * A section can inherit from another section by writing [section:
             * parentsection] - then all key-value-pairs of the parent section
             * are imported into section as default values
             */
            if (line[0] == '[') {
                size_t columnPos = line.find_first_of(':');
                if (columnPos == std::string::npos) {
                    sectionName = trim(line.substr(1, line.find_last_of(']') - 1));
                    if (sections.count(sectionName)) {
                        throw IniException("Multiple sections with the same name are not supported yet.");
                    }
                    sections[sectionName] = Description(sectionName);
                }
                else {
                    sectionName = trim(line.substr(1, columnPos-1));
                    parentSectionName = trim(line.substr(columnPos+1, line.find_last_of(']') - columnPos - 1));
                    if (sections.count(sectionName)) {
                        throw IniException("Multiple sections with the same name are not supported yet.");
                    }
                    if (!sections.count(parentSectionName)) {
                        throw IniException("Parent section '" + parentSectionName + "' does not exist.");
                    }
                    sections[sectionName] = Description(sectionName);
                    sections[sectionName].items = sections[parentSectionName].items;
                }

                continue;
            }
            /*
             * key = value pairs
             */
            size_t pos = line.find_first_of('=');
            if (pos == std::string::npos) {
                continue;
            }
            key = trim(line.substr(0, pos));
            value = trim(line.substr(pos+1));
            //skip empty keys, empty values are ok, though
            if (key == "") {
                continue;
            }
            //append it to the last section, if there is already a section
            if (sections.size() == 0) {
                throw IniException("Please specify at least one section.");
            }
            sections[sectionName].items[key] = DescriptionItem(value);
        }
    }

    std::string fileName;
};

class CommandLineOptionsException : public std::runtime_error
{
public:
    CommandLineOptionsException(const std::string& message)
        : std::runtime_error(message) {}
};

/*
 * TODO: keep order of command line arguments as they were specified (for
 * printing usage notice) 
 * TODO: properly format usage notice 
 * TODO: maybe make DescriptionItems available with both the long and the short
 * name, e.g. opts["hopping"] and ["t"]
 */
class CommandLineOptions: public Description
{
public:
    virtual ~CommandLineOptions() {}

    CommandLineOptions& add(Description option)
    {
        sections[option.name] = option;
        return *this;
    }

    CommandLineOptions& add(std::string name, std::string shortName, std::string description)
    {
        Description newOption(name);
        newOption["shortName"] = shortName;
        newOption["description"] = description;
        return add(newOption);
    }

    CommandLineOptions& add(std::string name, std::string description)
    {
        Description newOption(name);
        newOption["description"] = description;
        return add(newOption);
    }

    void parse(int argc, const char** argv)
    {
        for (int i=1; i<argc; i++) {
            std::string s(argv[i]);
            std::string optionName;
            //is it an option?
            if (s[0] == '-') {
                //is it a long or short option?
                if (s[1] == '-') {
                    optionName = s.substr(2);
                    if (!hasOption(optionName)) {
                        throw CommandLineOptionsException("Unknown option: '" + optionName + "'.");
                    }
                }
                else {
                    optionName = s.substr(1);
                    if (!hasShortOption(optionName)) {
                        throw CommandLineOptionsException("Unknown option: '" + optionName + "'.");
                    }
                    optionName = getOptionName(optionName);
                }
                //if there is a next argument and this argument isn't an option
                //specifier (i.e. not starting with '-') then this argument is the
                //value for our option
                if (i+1 < argc && argv[i+1][0] != '-') {
                    items[optionName] = DescriptionItem(std::string(argv[i+1]));
                    i++;
                }
                else {
                    items[optionName] = DescriptionItem("true");
                }
            }
            else {
                //we throw an error on unrecognised options / arguments
                throw CommandLineOptionsException("Can't parse command line arguments.");
            }
        }

    }

    std::string optionsDescription()
    {
        std::stringstream s;
        typedef SectionsType::iterator SIT;
        for (SIT i=sections.begin(); i!=sections.end(); i++) {
            if (i->second.hasItem("shortName")) {
                s << "-" << i->second.items["shortName"] << ", ";
            }
            s << "--" << i->second.name << "\t" << i->second.items["description"] << std::endl;
        }
        return s.str();
    }

private:
    bool hasOption(std::string name)
    {
        return hasSection(name);
    }

    bool hasShortOption(std::string shortName)
    {
        typedef SectionsType::iterator SIT;
        for (SIT i=sections.begin(); i!=sections.end(); i++) {
            if (i->second.hasItem("shortName") && i->second.items["shortName"] == shortName) {
                return true;
            }
        }
        return false;
    }

    std::string getOptionName(std::string shortName)
    {
        typedef SectionsType::iterator SIT;
        for (SIT i=sections.begin(); i!=sections.end(); i++) {
            if (i->second.hasItem("shortName") && i->second.items["shortName"] == shortName) {
                return i->second.name;
            }
        }
        return "";
    }
};

#endif /* __UTILITIES_HPP__ */
