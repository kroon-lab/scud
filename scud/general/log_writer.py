#!/usr/bin/env cctbx.python
import sys

class Log(object):

    def __init__(self,log_name='log.txt'):

        self.width = 80
        self.lf = open(log_name,'w')
        self.c_end = '\033[0m'
        
    def _print_hash_line(self):
        """
        Print 80 character long line with #
        """
        return '#'*self.width

    def _wrap_line(self,c=None,line=None):
        """
        Wrap ANSI characters around line
        """
        return c+line+self.c_end

    def close_log(self):
        self.lf.close()

    def process_message(self,line=None):
        """
        Write a process message
        """
        c = '\033[92m'
        proc_symbol = '**  '
        self.write(line = proc_symbol + line + '\n', c=c)

    def title(self,line=None):
        """
        Print line of text in between #'s
        """
        l = len(line)
        step = int((self.width - l - 2) / 2)
        s = '#' + ' '*step + line + ' '*step + '#'
        if len(s) > self.width: s[1]=''
        if len(s) < self.width: s = '# '+ s[1:]
        self.write(line='\n'+self._print_hash_line() + '\n' + s + '\n' + self._print_hash_line()+'\n\n')

    def show_info(self,line=None):
        """
        Write a info message
        """
        info_symbol = '  '
        self.write(info_symbol+line + '\n')

    def table(self, table_dict=None,length=15):
        """
        Print table of dictionairy

        dict should contain lists per element (same lenght)
        and one list names 'row' and one named 'column'
        """
        def _row(l):
            '''
            Format normal row
            '''
            s = ''
            for i in l : s += ' {i:{length}} |'.format(i=i,length=length)
            return s + '\n'
        def _2nd_row(l):
            '''
            Format non-text row
            '''
            s,ss = '',''
            for i in range(length): ss+='-'
            for i in range(l): s += ' {i:{length}} +'.format(i=ss,length=length)
            return s + '\n'
        
        # Get dict keys and remove row and column info
        key_list = table_dict.keys()
        key_list.remove('row')
        key_list.remove('column')

        f = ''
        f += _2nd_row(len(table_dict['column']))
        f += _row(table_dict['column'])
        f += _2nd_row(len(table_dict['column']))
        for row_num in range(len(table_dict[key_list[0]])):
            f += _row([table_dict['row'][row_num]] + [table_dict[k][row_num] for k in key_list])
        f += _2nd_row(len(table_dict['column']))
        self.write(line=f)
                
    def warning(self,line=None):
        """
        Write a warning message
        """
        c_warning = '\033[93m' # ANSI code
        warning_symbol = '!! '
        self.write(line=warning_symbol+line+'\n',c=c_warning)
        
    def write(self,line=None,c=None):
        """
        Write line to both screen and logfile
        """
        if c != None:
            sys.stdout.write(self._wrap_line(line=line,c=c))
        else:
            sys.stdout.write(line)
        sys.stdout.flush()
        self.lf.write(line)
        
        
        
