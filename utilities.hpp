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

#include <iostream>

#include <map>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <cctype>
#include <cstring>
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

class ConversionException : public std::runtime_error
{
public:
    explicit ConversionException (const std::string& message)
        : std::runtime_error(message) {}
};

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
        if (!s.eof() || value.empty())
            throw ConversionException("Can't convert '" + value + "' to double.");
        return returnValue;
    }

    operator int() const
    {
        int returnValue;
        std::stringstream s(value);
        s >> returnValue;
        if (!s.eof() || value.empty())
            throw ConversionException("Can't convert '" + value + "' to int.");
        return returnValue;
    }

    operator size_t() const
    {
        int returnValue = get<int>();
        if (returnValue < 0)
            throw ConversionException("Can't convert '" + value + "' to size_t.");
        return returnValue;
    }

    operator bool() const
    {
        std::string lvalue = toLower(value);
        if (lvalue == "no" || lvalue == "false")
            return false;
        else if (lvalue == "yes" || lvalue == "true")
            return true;

        throw ConversionException("Can't convert '" + value + "' to bool.");
    }

    operator std::string() const
    {
        return value;
    }
    
    /**
     * Access DescriptionItem as a list. 
     *
     * Convert DescriptionItem to a list (vector). Understands the formats
     * '1,2,3,4', '(1,2,3,4)', '1:2:0.1', '(1:2:0.1)' and also '1' (a list with
     * only one entry). Brackets are optional. A comma separated list is really
     * only a list, whereas three colon separated numbers constitute a
     * 'generator'. '1:2:0.1' generates a list starting at '1', ending at '2'
     * with '0.1' increments. Thus this example is equivalent to
     * '1.0,1.1,1.2,...,1.8,1.9,2.0'.
     *
     * @tparam T type of list
     *
     * @return list of type T
     */
    template<typename T>
    operator std::vector<T> ()
    {
        /*
         * TODO: obviously, std::vector won't scale up indefinitely; we could
         * use something similar to Python's generators instead
         */
        // brackets () are optional
        if (value.size() > 0 && value[0] == '(')
        {
            if (value.find(")") != value.size()-1)
                throw ConversionException("Missing closing bracket ')' in '" + value + "'.");
            value = value.substr(1, value.size()-2);
        }

        if (value.size() == 0)
            return std::vector<T>();

        if (value.find(',') != std::string::npos)
            return getList<T>(value, ',');
        
        if (value.find(':') != std::string::npos)
        {
            std::vector<T> tuple = getList<T>(value, ':');
            if (tuple.size() != 3) throw ConversionException("Invalid generator specification: '" + value + "'.");
            T begin = tuple[0];
            T end = tuple[1];
            T inc = tuple[2];
            T epsilon = inc * 0.0001;
            if (inc==0 || (begin<=end && inc<0) || (end<=begin && inc>0)) 
                throw ConversionException("Invalid generator specification: '" + value + "'.");
            std::vector<T> generated;
            for (T i=begin; (begin<=end && i<=end+epsilon) || 
                            (end<=begin && i>=end+epsilon); i+=inc)
                generated.push_back(i);
            return generated;
        }
        
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
    std::vector<T> getList (const std::string& v, const char sep) const
    {
        if (v.size() == 0) return std::vector<T>();
        std::vector<T> list;
        DescriptionItem item;
        size_t pos1=0;
        size_t pos2=0;
        while (pos1 <= v.size())
        {
            pos2 = v.find(sep, pos1);
            if (pos2 == std::string::npos) pos2 = v.size();
            if (pos2 == pos1) throw ConversionException("Error parsing list: '" + v + "'.");
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
                /*
                 * if there is a next argument and this argument isn't an
                 * option specifier (i.e. not a '-' followed by a non-digit)
                 * then this argument is the value for our option
                 */
                if (i+1 < argc &&   
                    ( (std::strlen(argv[i+1]) > 0 && argv[i+1][0] != '-') || 
                      (std::strlen(argv[i+1]) > 1 && std::isdigit(argv[i+1][1])) 
                    )
                   )
                {
                    items[optionName] = DescriptionItem(std::string(argv[i+1]));
                    i++;
                }
                else
                    items[optionName] = DescriptionItem("true");
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
