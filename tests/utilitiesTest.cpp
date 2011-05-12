/*
 * g++ -lboost_unit_test_framework -o utilitiesTest utilitiesTest.cpp
 */
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE utilities test
#include <boost/test/unit_test.hpp>

#include "utilities.hpp"
#include <cmath>

bool epsilonEqual (double a, double b)
{
    return std::fabs(a-b)<10E-10;
}

BOOST_AUTO_TEST_CASE ( empty_OptionValue )
{
    OptionValue i;
    BOOST_CHECK (i.getValue() == "");

    BOOST_CHECK_THROW (static_cast<double>(i), ConversionException);
    BOOST_CHECK_THROW (static_cast<int>(i), ConversionException);
    BOOST_CHECK_THROW (static_cast<size_t>(i), ConversionException);
    BOOST_CHECK (i == "");
    BOOST_CHECK_THROW (static_cast<bool>(i), ConversionException);
    BOOST_CHECK (i.isSet() == false);
    std::vector<int> v = i;
    BOOST_CHECK (v.size() == 0);
}

BOOST_AUTO_TEST_CASE ( OptionValue_type_conversion )
{
    OptionValue i("4.2");

    BOOST_CHECK (i == "4.2");
    BOOST_CHECK (static_cast<double>(i) == 4.2);
    BOOST_CHECK_THROW (static_cast<int>(i), ConversionException);
    BOOST_CHECK_THROW (static_cast<size_t>(i), ConversionException);
    std::vector<double> vd = i;
    BOOST_CHECK (vd.size() == 1);
    BOOST_CHECK (vd[0] == 4.2);
    BOOST_CHECK_THROW (std::vector<int> vi = i, ConversionException);

    i = -7.2;
    BOOST_CHECK (i == "-7.2");
    BOOST_CHECK (static_cast<double>(i) == -7.2);

    i = -5;
    BOOST_CHECK (i == "-5");
    BOOST_CHECK (static_cast<int>(i) == -5);
    BOOST_CHECK (static_cast<double>(i) == -5.0);
    BOOST_CHECK_THROW (static_cast<size_t>(i), ConversionException);

    i = 230;
    BOOST_CHECK (i == "230");
    BOOST_CHECK (static_cast<int>(i) == 230);
    BOOST_CHECK (static_cast<double>(i) == 230.0);
    BOOST_CHECK (static_cast<size_t>(i) == 230);

    i = "0";
    BOOST_CHECK_THROW (static_cast<bool>(i), ConversionException);

    i = "";
    BOOST_CHECK_THROW (static_cast<bool>(i), ConversionException);

    i = 4;
    BOOST_CHECK_THROW (static_cast<bool>(i), ConversionException);
    
    i = "false";
    BOOST_CHECK (static_cast<bool>(i) == false);
    BOOST_CHECK (static_cast<bool>(i) != true);

    i = "true";
    BOOST_CHECK (static_cast<bool>(i) == true);
    BOOST_CHECK (static_cast<bool>(i) != false);

    //TODO: test / check what happens when the number is too large for int or
    //size_t, etc
}

BOOST_AUTO_TEST_CASE ( OptionValue_default_value_and_isSet )
{
    OptionValue i;

    BOOST_CHECK (i.isSet() == false);
    BOOST_CHECK_THROW (static_cast<bool>(i), ConversionException);

    i.setDefault(7);
    BOOST_CHECK (i.isSet() == false);
    BOOST_CHECK (static_cast<int>(i) == 7);

    i.setDefault(true);
    BOOST_CHECK (i.isSet() == false);
    BOOST_CHECK (static_cast<bool>(i) == true);

    i.setDefault(-8.4);
    BOOST_CHECK (i.get(5) == 5);
    BOOST_CHECK (i.get<double>() == -8.4);
    BOOST_CHECK (i.isSet() == false);

    i = -8.4;
    BOOST_CHECK (i.get<double>(5) == -8.4);
    BOOST_CHECK (i.get<double>() == -8.4);
    BOOST_CHECK (i.isSet() == true);
}

BOOST_AUTO_TEST_CASE ( OptionValue_lists )
{
    OptionValue i;
    std::vector<double> v;

    i = "[]";
    BOOST_CHECK_THROW (v = i, ConversionException);
    BOOST_CHECK (v.size() == 0);

    i = "[10,20,30,40,50]";
    BOOST_CHECK_THROW (v = i, ConversionException);
    
    i = "";
    v = i;
    BOOST_CHECK (v.size() == 0);

    i = "()";
    v = i;
    BOOST_CHECK (v.size() == 0);

    i = "(1)";
    v = i;
    BOOST_CHECK (v.size() == 1);
    BOOST_CHECK (v[0] == 1);

    i = "1";
    v = i;
    BOOST_CHECK (v.size() == 1);
    BOOST_CHECK (v[0] == 1);

    i = "10,20,30,40,50";
    v = i;
    BOOST_CHECK (v.size() == 5);
    BOOST_CHECK (v[3] == 40);

    i = "(10,20,30,40,50,)";
    BOOST_CHECK_THROW (v = i, ConversionException);
    i = "(10,20,30,40,50";
    BOOST_CHECK_THROW (v = i, ConversionException);
    i = "10,20,30,40,50)";
    BOOST_CHECK_THROW (v = i, ConversionException);
    i = "10,20,(30,40,50)";
    BOOST_CHECK_THROW (v = i, ConversionException);
    i = "(10,20),30,40,50";
    BOOST_CHECK_THROW (v = i, ConversionException);
    i = "10,20,30,40,50,";
    BOOST_CHECK_THROW (v = i, ConversionException);
    i = "10,20,30,40,,50";
    BOOST_CHECK_THROW (v = i, ConversionException);
    i = "10,20,30,4 0,50";
    BOOST_CHECK_THROW (v = i, ConversionException);

    i = "10, 20,30,40 ,50";
    v = i;
    BOOST_CHECK (v.size() == 5);
    BOOST_CHECK (v[3] == 40);

    i = "10:20";
    BOOST_CHECK_THROW (v = i, ConversionException);
    i = "10:20:30:40";
    BOOST_CHECK_THROW (v = i, ConversionException);

    i = "2:1:0.1";
    BOOST_CHECK_THROW (v = i, ConversionException);
    i = "1:2:-0.1";
    BOOST_CHECK_THROW (v = i, ConversionException);
    i = "(1:2:-0.1)";
    BOOST_CHECK_THROW (v = i, ConversionException);
    i = "2:1:0.1:";
    BOOST_CHECK_THROW (v = i, ConversionException);
    
    i = "(1:2:0.1)";
    v = i;
    BOOST_CHECK (v.size() == 11);
    BOOST_CHECK (epsilonEqual(v[2], 1.2));
    
    i = "1:2:0.1";
    v = i;
    BOOST_CHECK (v.size() == 11);
    BOOST_CHECK (epsilonEqual(v[2], 1.2));
    
    i = "2:1:-0.1";
    v = i;
    BOOST_CHECK (v.size() == 11);
    BOOST_CHECK (epsilonEqual(v[2], 1.8));
}

BOOST_AUTO_TEST_CASE ( OptionSection_addSection )
{
    OptionSection os("Main");
    os.addSection("Part1");
    os.addSection("Part2");
    BOOST_CHECK_THROW (os.addSection("Part1"), OptionSectionException);

    BOOST_CHECK (os.hasSection("Part1") == true);
    BOOST_CHECK (os.hasSection("Part2") == true);
    BOOST_CHECK (os.hasSection("Part3") == false);
}

BOOST_AUTO_TEST_CASE ( OptionSection_addOption )
{
    OptionSection os("Main");
    os.addSection("Part1");

    os.s("Part1").o("a") = 0.1;
    os.s("Part1")
        .addOption("b")
        .addOption("c");
    
    BOOST_CHECK_THROW (os.s("Part1").addOption("a"), OptionSectionException);
    BOOST_CHECK (os.s("Part1").o("b") == "");
    double a = os.s("Part1")["a"];
    BOOST_CHECK (epsilonEqual(a, 0.1));
    a = os.getSection("Part1").getOption("a").getValue();
    BOOST_CHECK (epsilonEqual(a, 0.1));
}

BOOST_AUTO_TEST_CASE ( OptionSection_nestedSections )
{
    OptionSection os("Main");
    os.addSection("Part1");

    os["a"] = 5;
    os["b"] = "huhu";
    
    os.s("Part1")["a"] = 1.7;
    os.s("Part1").o("b") = "blub";

    os.s("Part1").addSection("Part1A");
    os.s("Part1").s("Part1A")["a"] = "blip";
    os.s("Part1").s("Part1A")["b"] = 10000;

    BOOST_CHECK (static_cast<int>(os["a"]) == 5);
    BOOST_CHECK (os.o("b") == "huhu");
    BOOST_CHECK (static_cast<double>(os.s("Part1")["a"]) == 1.7);
    BOOST_CHECK (os.s("Part1")["b"] == "blub");
    BOOST_CHECK (os.s("Part1").s("Part1A").o("a") == "blip");
    BOOST_CHECK (static_cast<int>(os.s("Part1").s("Part1A").getOption("b").getValue()) == 10000);
}

CommandLineOptions setupCLOptions ()
{
    CommandLineOptions o;
    o.add("help", "h", "Print this help message.")
     .add("model", "m", "Which QCA model to use.")
     .add("N_p", "p", "Number of plaquets.")
     .add("", "t", "Hopping parameter.") // only short option, i.e. -t, but not --t
     .add("", "a", "Intra-plaquet spacing.")
     .add("", "b", "Inter-plaquet spacing.")
     .add("energy-spectrum", "E", "Calculate the energy spectrum.")
     .add("polarisation", "P", "Calculate the polarisation for specified plaquet(s).");
    
    o["N_p"].setDefault(1);
    o["t"].setDefault(1);
    o["b"].setDefault(1.75);

    return o;
}

BOOST_AUTO_TEST_CASE ( CommandLineOptions_setupCLO )
{
    CommandLineOptions o = setupCLOptions();

    BOOST_CHECK (o.hasOption("model") == false); //because it has no default value
    BOOST_CHECK (o.hasOption("N_p") == true); //because it has a default value
    BOOST_CHECK (o.hasOption("humpa") == false); //because this option does not exist
    BOOST_CHECK (static_cast<int>(o["t"]) == 1);
    BOOST_CHECK (epsilonEqual(o["b"], 1.75));
}

BOOST_AUTO_TEST_CASE ( CommandLineOptions_long_and_short_options )
{
    CommandLineOptions o = setupCLOptions();

    o["model"] = "Blub";
    o["polarisation"] = true;
    o["energy-spectrum"] = false;
    o["t"] = 1;

    BOOST_CHECK (o["m"] == "Blub");
    BOOST_CHECK (static_cast<bool>(o["P"]) == true);
    BOOST_CHECK (static_cast<bool>(o["E"]) == false);
    BOOST_CHECK (epsilonEqual(o["t"], 1));

    BOOST_CHECK (o.hasOption("model") == true);
    BOOST_CHECK (o.hasOption("m") == true);
    BOOST_CHECK (o.hasOption("polarisation") == true);
    BOOST_CHECK (o.hasOption("P") == true);
    BOOST_CHECK (o.hasOption("t") == true);
    BOOST_CHECK (o.hasOption("a") == false);
}

void fillArgv (char** argv, const char* s0="", const char* s1="", const char* s2="", 
                const char* s3="", const char* s4="", const char* s5="", 
                const char* s6="", const char* s7="", const char* s8="", 
                const char* s9="")
{
    strcpy(argv[0], s0);
    strcpy(argv[1], s1);
    strcpy(argv[2], s2);
    strcpy(argv[3], s3);
    strcpy(argv[4], s4);
    strcpy(argv[5], s5);
    strcpy(argv[6], s6);
    strcpy(argv[7], s7);
    strcpy(argv[8], s8);
    strcpy(argv[9], s9);
}

//void fillArgv (char** argv, const char* s0="", const char* s1="", const char* s2="", 
//                const char* s3="", const char* s4="", const char* s5="", 
//                const char* s6="", const char* s7="", const char* s8="", 
//                const char* s9="")
//{
//    strcpy(*argv, s0);
//    argv += std::strlen(s0) + 1;
//    strcpy(*argv, s1);
//    argv += std::strlen(s0) + 1;
//    strcpy(*argv, s2);
//    argv += std::strlen(s0) + 1;
//    strcpy(*argv, s3);
//    argv += std::strlen(s0) + 1;
//    strcpy(*argv, s4);
//    argv += std::strlen(s0) + 1;
//    strcpy(*argv, s5);
//    argv += std::strlen(s0) + 1;
//    strcpy(*argv, s6);
//    argv += std::strlen(s0) + 1;
//    strcpy(*argv, s7);
//    argv += std::strlen(s0) + 1;
//    strcpy(*argv, s8);
//    argv += std::strlen(s0) + 1;
//    strcpy(*argv, s9);
//}

BOOST_AUTO_TEST_CASE ( CommandLineOptions_parser )
{
    int argc;
    char** argv = new char* [10];
    for (size_t i=0; i<10; i++)
        argv[i] = new char [100];

    CommandLineOptions o = setupCLOptions();
    fillArgv(argv, "", "--m", "blub");
    argc = 3;
    BOOST_CHECK_THROW (o.parse(argc, const_cast<const char**>(argv)), CommandLineOptionsException);

    o = setupCLOptions();
    fillArgv(argv, "", "--t", "4.3");
    argc = 3;
    BOOST_CHECK_THROW (o.parse(argc, const_cast<const char**>(argv)), CommandLineOptionsException);

    o = setupCLOptions();
    fillArgv(argv, "", "--a", "2");
    argc = 3;
    BOOST_CHECK_THROW (o.parse(argc, const_cast<const char**>(argv)), CommandLineOptionsException);

    o = setupCLOptions();
    fillArgv(argv, "", "-N_p", "1");
    argc = 3;
    BOOST_CHECK_THROW (o.parse(argc, const_cast<const char**>(argv)), CommandLineOptionsException);

    o = setupCLOptions();
    fillArgv(argv, "", "--bogus");
    argc = 2;
    BOOST_CHECK_THROW (o.parse(argc, const_cast<const char**>(argv)), CommandLineOptionsException);

    o = setupCLOptions();
    fillArgv(argv, "", "-bogus");
    argc = 2;
    BOOST_CHECK_THROW (o.parse(argc, const_cast<const char**>(argv)), CommandLineOptionsException);

    o = setupCLOptions();
    fillArgv(argv, "", "-m", "blub", "-t", "4.3", "-P", "--N_p", "1");
    argc = 8;
    o.parse(argc, const_cast<const char**>(argv));
    BOOST_CHECK (o.hasOption("model") == true);
    BOOST_CHECK (o["model"] == "blub");
    BOOST_CHECK (o.hasOption("t") == true);
    BOOST_CHECK (epsilonEqual(o["t"], 4.3));
    BOOST_CHECK (o.hasOption("polarisation") == true);
    BOOST_CHECK (static_cast<bool>(o["polarisation"]) == true);
    BOOST_CHECK (static_cast<int>(o["p"]) == 1);
    BOOST_CHECK (o.hasOption("a") == false);
    BOOST_CHECK (o.hasOption("b") == true);
    BOOST_CHECK (epsilonEqual(o["b"], 1.75));

    o = setupCLOptions();
    fillArgv(argv, "", "--model", "blub", "-t", "4.3", "--polarisation", "-p", "1");
    argc = 8;
    o.parse(argc, const_cast<const char**>(argv));
    BOOST_CHECK (o.hasOption("model") == true);
    BOOST_CHECK (o["model"] == "blub");
    BOOST_CHECK (o.hasOption("t") == true);
    BOOST_CHECK (epsilonEqual(o["t"], 4.3));
    BOOST_CHECK (o.hasOption("polarisation") == true);
    BOOST_CHECK (static_cast<bool>(o["polarisation"]) == true);
    BOOST_CHECK (static_cast<int>(o["p"]) == 1);
    BOOST_CHECK (o.hasOption("a") == false);
    BOOST_CHECK (o.hasOption("b") == true);
    BOOST_CHECK (epsilonEqual(o["b"], 1.75));

    //TODO: how do we handle duplicate command line options? i.e. -t 2 -t 4
}

BOOST_AUTO_TEST_CASE ( CommandLineOptions_usage_message )
{
}
