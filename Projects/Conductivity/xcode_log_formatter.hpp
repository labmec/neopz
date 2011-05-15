#ifndef XCODE_LOG_FORMATTER_HPP_INCLUDED
#define XCODE_LOG_FORMATTER_HPP_INCLUDED

#include <boost/test/output/compiler_log_formatter.hpp>

struct xcode_log_formatter : 
	public boost::unit_test::output::compiler_log_formatter
{
	void print_prefix(std::ostream& output, 
		boost::unit_test::const_string file_name, std::size_t line)
	{
		output << file_name << ':' << line << ": ";
	}
};

#endif // XCODE_LOG_FORMATTER_HPP_INCLUDED
