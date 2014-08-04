/*
 * minimal code
 *
 
 #include "CommandLineParser.hpp"
    p::CommandLineParser parser(argc, argv);

    int input = parser.addOption<int>("-i",-17,"test int");
    std::string allo = parser.addOption<std::string>("-s","coucou","test string");
    std::string peep = parser.addOption<std::string>("-p","salut");
 
 *
 * To do
 *
  
  add default -h option
 
 */
#include <iostream>
#include <cstdlib> //for EXIT_SUCCESS
#include <sstream>
#include <string>
#include <vector>
#include <map>

namespace p
{
    // STRING INT
    void ConvertType(std::string& s, int& i)
	{
		i = std::atoi( s.c_str() );
	}
    
    void ConvertType(int& i, std::string& s)
    {
        std::stringstream ss;
        ss << i;
        ss >> s;
    }
    
    // STRING FLOAT
	void ConvertType(std::string& s, float& f)
	{
		f = (float)std::atof( s.c_str() );
	}
    
    void ConvertType(float& f, std::string& s)
    {
        std::stringstream ss;
        ss << f;
        ss >> s;
    }
    
    // STRING DOUBLE
	void ConvertType(std::string& s, double& d)
	{
		d = std::atof( s.c_str() );
	}
    
    void ConvertType(double& d, std::string& s)
    {
        std::stringstream ss;
        ss << d;
        ss >> s;
    }
	
    // STRING CHAR*
	void ConvertType(std::string& s, char* c)
	{
        //c_str returns a const char*
		c = const_cast<char*>( s.c_str() );
	}
    
    void ConvertType(char* c, std::string& s)
    {
        std::stringstream ss;
        ss << c;
        ss >> s;
    }
    
    
    
class CommandLineParser
{
    
private:
    
int								   _numArg;
std::vector<std::string>		   _argList;
std::map<std::string, std::string> _argmap;
std::map<std::string, std::string> _descriptionList;


public:
CommandLineParser(int argc, char** argv)
{
	_numArg = argc;

	// convert char** to vector<string>

	for (int i = 1; i < argc; i++)
	{
		std::stringstream ss;
		std::string		  s;
		ss << argv[i];
		ss >> s;
		_argList.push_back(s);
        
        if ( s =="-h" || s =="--help")
        {
            HelpMessage();
            exit(EXIT_SUCCESS);
        }
        
        
	}

	// map <option , Argument>

	for (int i = 0; i < _argList.size(); i++)
	{
		if (_argList[i][0] == '-')
		{
			if (_argList[i + 1][0] != '-')
				_argmap.insert(std::pair<std::string, std::string>(_argList[i], _argList[i + 1]) );
			else
				_argmap.insert(std::pair<std::string, std::string>(_argList[i], "") );
		}
	}
}

template<typename Type>
Type addOption(std::string optName, Type defaultValue, std::string description = "default description")
{
	Type result;

	//if option is found

	if (_argmap.find(optName) != _argmap.end() )
    {
		//ConvertType(_argmap.find(optName)->second, result);
        std::stringstream ss;
        ss << _argmap.find(optName)->second;
        ss >> result;
        
    }
	//if option not found, assign default value

	else
		result = defaultValue;
    
    //add entry to descriptionList
    _descriptionList.insert(std::pair<std::string, std::string>(optName, description) );
    

	return result;
}
    
    void HelpMessage()
    {
/*        std::map<std::string, std::string>::iterator it_descr, it_arg;
        
        std::cout<<"Usage :"<<std::endl;
        for (it_descr=_descriptionList.begin(), it_arg=_argmap.begin();
             it_descr<_descriptionList.end();
             it_descr++, it_arg++)
        {
            std::cout<<"\t"<<it_descr->first<<" :\t"<<it_descr->second<<"\t[default = "<<it_arg->find(it_descr->first)->second<<"]"<<std::endl;
        }
 */
    }
    
    
};


/*
 *
 * CommandLineParser parser(argc, argv);
 * parser.addOption<int>( "-i", &input);
 * parser.addOption<std::string>( "-s", &allo, &quiEst, &la );
 *
 *
 *** overload
 *
 *
 * int i = parser.addOption<int>("-i", defaultValue);
 * std::string s[] = parser.addOption<std::string>( "-s", "defaultValue" );
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * */
}
