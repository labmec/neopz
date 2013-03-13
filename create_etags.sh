find ./ -regex ".*[cChH]\(pp\)?" -print  | grep -v svn | grep -v OutOfDate | etags --language=c++ -
