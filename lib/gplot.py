# A python-gnuplot interface

from __future__ import print_function
import subprocess
from numpy import savetxt
import tempfile
import os

__author__ = 'Mathias Zechmeister'
__version__ = 'v13'
__date__ = '2018-10-25'
__all__ = ['gplot', 'Gplot', 'ogplot', 'Iplot']


class Gplot(object):
   """
   An interface between Python and gnuplot.
   
   Creation of an instance opens a pipe to gnuplot and return an object for communication.
   Plot commands are send to gnuplot via the call method; arrays as arguments are handled .
   Gnuplot options are set by calling them as method attributes.
   Each method returns the object again. This allows to chain set and plot method.
   
   Parameters
   ----------
   tmp : str, optional
       Method for passing data.
       * None - create a non-persistent temporary file (default)
       * '' - create a local persistent file
       * '-' - use gnuplot special filename (no interactive zoom available,
               replot does not work)
       * '$' - use inline datablock (not faster than temporary data,
             does not work with flush='' an ogplot)
       * 'filename' - create manually a temporary file
   stdout : boolean, optional
       If true, plot commands are send to stdout instead to gnuplot pipe.
   mode : str, optional
       Primary command for the call method. The default is 'plot'. After creation it can
       be changed, e.g. gplot.mode = gplot.splot.

   args : array or str for function, file, or other plot commands like style
   flush : str, optional
       set to '' to suppress flush until next the ogplot (for large data sets)

   Methods
   -------
   __call__
   load
   replot
   plot
   print
   put
   set
   splot
   test
   unset
   var
   
   NOTES
   -----
   The attribute print does not work in python 2.

   Examples
   --------
   
   A simple plot and add a data set
   
   >>> gplot('sin(x) w lp lt 2')
   >>> gplot+(np.arange(10)**2., 'w lp lt 3')
   >>> gplot+('"filename" w lp lt 3')
   
   Pass multiple curves in one call
   
   >>> gplot('x/5, 1, x**2/50 w l lt 3,', np.sqrt(np.arange(10)),' us 0:1 ps 2 pt 7 ,sin(x)')
   >>> gplot.mxtics().mytics(2).repl
   >>> gplot([1,2,3,4])
   >>> gplot([1,2,3,4], [2,3,1,1.5])
   >>> gplot([1,2,3,4], [[2,2,1,1.5], [3,1,4,5.5]])
   >>> gplot([[2,2,1,1.5]])
   >>> gplot([1],[2],[3],[4])
   >>> gplot(1,2,3,4)
   
   """
   version = subprocess.check_output(['gnuplot', '-V'])
   version = float(version.split()[1]) 
   
   def __init__(self, cmdargs='', tmp=None, mode='plot', stdout=False):
      self.stdout = stdout
      self.tmp = tmp
      self.mode = getattr(self, mode)   # set the default mode for __call__ (plot, splot)
      self.gnuplot = subprocess.Popen('gnuplot '+cmdargs, shell=True, stdin=subprocess.PIPE,
                   universal_newlines=True, bufsize=0)   # This line is needed for python3! Unbuffered and to pass str instead of bytes
      self.pid = self.gnuplot.pid
      self.og = 0   # overplot number
      self.buf = ''
      self.tmp2 = []
      self.flush = None
   
   def _plot(self, *args, **kwargs):
      # collect all arguments
      tmp = kwargs.pop('tmp', self.tmp)
      flush = kwargs.pop('flush', '\n')
      if self.version in [4.6] and flush=="\n": flush = "\n\n"   # append a newline to workaround a gnuplot pipe bug
      # with mouse zooming (see http://sourceforge.net/p/gnuplot/bugs/1203/)
      self.flush = flush
      pl = ''
      data = ()
      for arg in args + (flush,):
         if isinstance(arg, str):   # append argument, but flush the data before
            if data:
               # transpose data when writing
               data = zip(*data)
               self.og += 1
               tmpname = tmp
               if tmp in ('-',):
                  # use gnuplot's special filename '-'
                  self.buf += "\n".join(" ".join(map(str,tup)) for tup in data)+"\ne\n"
               elif tmp in ('$',):
                  # gnuplot's inline datablock
                  tmpname = "$data%s" % self.og
                  # prepend the datablock
                  pl = tmpname+" <<EOD\n"+("\n".join(" ".join(map(str,tup)) for tup in data))+"\nEOD\n" + pl
               elif tmp is None:
                  # create temporary file; default
                  self.tmp2.append(tempfile.NamedTemporaryFile())
                  tmpname = self.tmp2[-1].name
                  #savetxt(tmp2[-1], data) #zipped data not supported by python3
                  savetxt(self.tmp2[-1], list(data), fmt="%s") # zipped data not supported by python3, fmt="%s" to allow for strings
                  self.tmp2[-1].seek(0)
               else:
                  # create local temporary file
                  if tmp=='':
                     tmpname = 'gptmp_'+str(self.pid)+str(self.og)
                  savetxt(tmpname, data, fmt="%s")
               pl += '"'+tmpname+'"'
            pl += arg
            data = ()
         else:
            # collect data; append columns and matrices
            _1D = hasattr(arg, '__iter__')
            _2D = _1D and hasattr(arg[0], '__iter__') and not isinstance(arg[0], str)
            data += tuple(arg) if _2D else (arg,) if _1D else ([arg],)
      self.put(pl, end='')
      if flush!='': self.put(self.buf, end='')

   def put(self, *args, **kwargs):
      # send the commands to gnuplot
      print(file=None if self.stdout else self.gnuplot.stdin, *args, **kwargs)
      return self

   # some plot commands (kwargs possible)
   def plot(self, *args, **kwargs):
      self.og = 0; self.buf = ''; self.put('\n')        # reset
      return self._plot('plot ', *args, **kwargs)

   def splot(self, *args, **kwargs):
      self.og = 0; self.buf = ''; self.put('\n')        # reset
      return self._plot('splot ', *args, **kwargs)

   def replot(self, *args, **kwargs):
      return self._plot('replot ', *args, **kwargs)

   def test(self, *args, **kwargs):
      return self._plot('test', *args, **kwargs)

   def oplot(self, *args, **kwargs):
      pl = ',' if self.flush=='' else ' replot '
      return self._plot(pl, *args, **kwargs)

   def var(self, **kwargs):
      # set gnuplot variables
      for i in kwargs.items(): self.put("%s=%s" % i)
      return self

   def __call__(self, *args, **kwargs):
      # by default plot mode is executed, but the user can change that
      return self.mode(*args, **kwargs)

   def __getattr__(self, name):
      # generic translatation, e.g. gplot.title sends "set title"
      if name in ('__repr__', '__str__'):
         raise AttributeError
      elif name=='repl':
         return self.replot()
      elif name in ['load', 'set', 'show', 'unset', 'reset', 'print', 'bind']:
         # some fixed keywords
         # print as attribute does not work in python 2
         def func(*args):
            return self.put(name, *args)
         return func
      else:
         # dynamic attributes (xlabel, key, etc.)
         def func(*args):
            return self.set(name, *args)
         return func

   def __add__(self, other):
      '''Add a curve to an existing plot (similar as ogplot).
      
      This will look like calling the gplot instance, but a tuple is retrieved.
      Hence no keywords can be passed.
      The plot will be update immediately. To pass flush='', use the __lt__
      method.
      '''
      self.oplot(*other)

   def __sub__(self, other):
      '''Start a new plot, but do not flush.
      
      In particular for larger data sets, one may want to accumulate
      remaing curves before, instead of replotting each time.
      Use gplot< to append more sets and gplot+ to finish the plot.
      The "-" may remind to clear previous plots or start a line of new curves.
      '''
      self(*other, flush='')

   def __lt__(self, other):
      '''Add a curve to an existing plot, but do not flush.
      
      In the "<" sign shall remind to the bash pipe operator.
      '''
      self.oplot(*other, flush='')


class Iplot(Gplot):
   '''
   Gnuplot for Jupyter.
   
   The class is similar as Gplot, but, instead of self, plot method returns
   an object that can be displayed.
   
   '''
   # Search gnuplot installation with javascript library
   # '/usr/share/gnuplot/%s/js/' % Gplot.version
   # '/usr/local/share/gnuplot/%s/js/' % Gplot.version
   # The path can be extracted from output of term svg mousing.
   # Request a small plot and parse the output
   _jsdir = subprocess.check_output(['gnuplot', '-e', 'set term svg mousing; unset border; unset tics; set samp 2; pl [][0:1]NaN t ""'], universal_newlines=True).split('gnuplot_svg.js')[0].split('xlink:href="')[-1]
   # universal_newlines=True leads finally to the same type|class=str in python 2 and python 3, instead of bytes

   # The absolute path might not accessible in Jupyter, if localhost was
   # started in local directory.

   # online locations:
   # latest version, including mouse wheel for zooming for svg:
   _jsdir = "http://gnuplot.sourceforge.net/demo_canvas_cvs/" 
   # "http://gnuplot.sourceforge.net/demo_canvas_5.2/"
   # "http://gnuplot.sourceforge.net/demo_svg/"
   
   def __init__(self, *args, **kwargs):
      '''
      Iplot(*args, suffix='png', uri=True, cleanup=True, **kwargs)
      
      Parameters:
      -----------
      *args : cmdargs (see Gplot).
      suffix : str, optional
          Default is 'png'.
          Support for svg, html, and png.
          svg looks nice. canvas (js) is more interactive.
      opt : str, optional
          Additional terminal settings (e.g. 'size 300,200').
      uri : boolean, optional
          If true, the figure will inline. A fifo is used, leading to a blocking
          read and has no problem with waiting for termination of an asynchronous
          process. It seems to be faster.
      cleanup : boolean, optional
          If true (default), the temporary image file will be deleted. For uri=False,
          this may lead to non-existing file.
      jsdir : str, optional
          Path to gnuplot javascript library.
          (default: http://gnuplot.sourceforge.net/demo_canvas_cvs/)
      **kwargs : tmp (see Gplot).
      
      Examples
      --------
      >>> iplot = Iplot(opt='size 300,200')
      >>> iplot('x')
      
      NOTES
      -----
      To get nice working gnuplot output, some issues needed to be
      solved.
         * canvas: To enable correct mousing, div#site position (which defaults
             to static) is set to sticky. Then it becomes an offsetParent and
             its scrollTop value resulting from overflow:auto can be processed
             in gnuplot_mouse.js, which was modified.
         * svg: Inline javascript is not loaded when using uri.
         * png: ok.
      Changing size was not tested yet.
      
      '''
      self.suffix = kwargs.pop('suffix', 'png')
      self.opt = kwargs.pop('opt', '')
      self.uri = kwargs.pop('uri', True)
      self.jsdir = kwargs.pop('jsdir', 'jsdir "%s"'% self._jsdir)
      self.cleanup = kwargs.pop('cleanup', True)
      self.canvasnum = 0
      return super(Iplot, self).__init__(*args, **kwargs)

   def _plot(self, *args, **kwargs):
      self.canvasnum += 1
      uri = self.uri
      canvasname = "fishplot_%s" % self.canvasnum
      term = {'svg': 'svg mouse %s' % self.jsdir,
              'svg5': 'svg mouse standalone name "bla"',
              'html': 'canvas name "%s" mousing %s' % (canvasname, self.jsdir)
             }.get(self.suffix, 'pngcairo') + ' ' + self.opt
      # png + no_uri, Image still converts into uri -> but works
      # png + uri, Image still converts into uri -> works
      # svg + local_svg, mousing works
      # svg + no_uri, as inline
      # svg + uri, as inline
      # html + no_uri (if cleanup=False, the browser can read the file)
      # html + uri -> works
      imgfile = tempfile.NamedTemporaryFile(suffix='.'+self.suffix).name
      if self.suffix=='svg' and not uri:
         # use a local file
         imgfile = 'simple_%s.svg' % self.canvasnum
      if self.suffix=='html':
         imgfile = canvasname + '.js'

      if uri:
         # cleanup if already exists
         os.system("rm -f "+imgfile)
         os.mkfifo(imgfile)

      self.term(term).out('"%s"' % imgfile)
      if uri:
         # the fifo needs something to read; but display will finally open and read imgfile
         fifo = open(imgfile, 'r')
      super(Iplot, self)._plot(*args, **kwargs)
      self.out()

      if uri:
         1
         #imgdata = fifo.read()
         #imgfile = 'data:image/png;base64,'+
         #imgdata = open(imgfile, "rb").read() # png needs rb, but imgfile can be also passed directly 
         #imgdata.replace('<svg ','<script type="text/javascript" xmlns:xlink="http://www.w3.org/1999/xlink" xlink:href="./gp/gnuplot_svg.js"></script>\n<svg ')
         #import base64
         #print (base64.b64encode(imgdata).decode('ascii')[-10:])
      if not uri:
         import time
         counter = 0
         while counter<100 and not (os.path.exists(imgfile) and os.system("lsof "+imgfile)):
            time.sleep(0.003)
            counter += 1
            if not self.cleanup: print(counter, end='\r')
         if not self.cleanup:
            print(counter, imgfile, os.path.exists(imgfile), os.system("lsof "+imgfile))

      imgdata = imgfile

      from IPython.display import Image, SVG, HTML, Javascript
      showfunc = {'svg':SVG, 'html':HTML}.get(self.suffix, Image)

      if self.suffix=='svg':
         if not uri:
            showfunc = HTML
            imgdata = '<embed src="%s" type="image/svg+xml">' % imgfile
         else:
            showfunc = HTML
            imgdata = open(imgfile).read() # .replace("onload","onclick")
            # to get mousing in chrome uncomment in gnuplot_svg.js:63
            # // p.x = evt.pageX; p.y = evt.pageY;
            if self._jsdir.startswith("http://"):
               from urllib.request import urlopen
               gnuplot_svg_js = urlopen(self._jsdir+"gnuplot_svg.js").read().decode('utf-8')
            else:
               gnuplot_svg_js = open(self._jsdir+"gnuplot_svg.js").read()
            gnuplot_svg_js.replace("documentElement", 'getElementsByTagName("svg")[0]')
            imgdata = ('''
    Inline SVG1
     <script>
''' + gnuplot_svg_js + '''
        // manually gnuplot_svg.Init (onload does not fire?)
        gnuplot_svg.SVGDoc = document.getElementsByTagName("svg")[0];
        console.log(gnuplot_svg.SVGRoot, gnuplot_svg.SVGDoc);
     </script>
    <object>
            ''' + imgdata + '</object>')

      if self.suffix=='html':
         if uri:
            imgdata = '<script>%s</script>' % open(imgfile).read()
         else:
            imgdata = '''<script src="%s"></script>''' % imgfile
         # style the buttons and work around some mousing issues 
         imgdata = '''
         <link rel="stylesheet" href="%sgnuplot_mouse.css" type="text/css">
<style>
   img.icon-image {
      max-width: None !important;
   }
   td.icon, td.mb0, td.mb1 {
      background-color: #f7f7f7;
      line-height: 16px;
   }
</style>
<script src="%scanvastext.js"></script>
<script src="%sgnuplot_common.js"></script>
<script src="%sgnuplot_mouse.js"></script>
<script type="text/javascript">
   // to get correct mouse coords (overflow:auto, scrollTop)
   document.getElementById("site").style.position = "sticky";
</script>
<script type="text/javascript">
   // from view-source:http://gnuplot.sourceforge.net/demo_canvas/simple.html
   var canvas, ctx;
   gnuplot.grid_lines = true;
   gnuplot.zoomed = false;
   gnuplot.active_plot_name = "gnuplot_canvas";
   gnuplot.active_plot = gnuplot.dummyplot;
   gnuplot.help_URL = "http://gnuplot.sourceforge.net/demo_canvas_5.0/canvas_help.html";
   gnuplot.dummyplot = function() {};
   function gnuplot_canvas( plot ) { gnuplot.active_plot(); };
</script>
         %s
<table class="noborder" style="margin-top:0; border:0;">
  <tr style="border:0;">
    <td style="border:0;">
    <table class="mbright" tabindex=0>
    <tr>
      <td class="icon" onclick="gnuplot.toggle_grid();"><img src="%sgrid.png" id="gnuplot_grid_icon" class="icon-image" alt="#" title="toggle grid"></td>
      <td class="icon" onclick="gnuplot.unzoom();"><img src="%spreviouszoom.png" id="gnuplot_unzoom_icon" class="icon-image" alt="unzoom" title="unzoom"></td>
      <td class="icon" onclick="gnuplot.rezoom();"><img src="%snextzoom.png" id="gnuplot_rezoom_icon" class="icon-image" alt="rezoom" title="rezoom"></td>
      <td class="icon" onclick="gnuplot.toggle_zoom_text();"><img src="%stextzoom.png" id="gnuplot_textzoom_icon" class="icon-image" alt="zoom text" title="zoom text with plot"></td>
      <td class="icon" onclick="gnuplot.popup_help();"><img src="%shelp.png" id="gnuplot_help_icon" class="icon-image" alt="?" title="help"></td>
      <td class="icon"></td>
      <td class="mb0">x&nbsp;</td> <td class="mb1"><span id="%s_x">&nbsp;NaN&nbsp;</span></td>
      <td class="mb0">y&nbsp;</td> <td class="mb1"><span id="%s_y">&nbsp;NaN&nbsp;</span></td>
      <td class="icon">&nbsp;</td>
      <td class="icon" onclick='gnuplot.toggle_plot("%s_plot_1")'>&#9312;</td>
      <td class="icon" onclick='gnuplot.toggle_plot("%s_plot_2")'>&#9313;</td>
      <td class="icon" onclick='gnuplot.toggle_plot("%s_plot_3")'>&#9314;</td>
      <td class="icon" onclick='gnuplot.toggle_plot("%s_plot_4")'>&#9315;</td>
      <td class="icon" onclick='gnuplot.toggle_plot("%s_plot_5")'>&#9316;</td>
      <td class="icon" onclick='gnuplot.toggle_plot("%s_plot_6")'>&#9317;</td>
      <td class="icon" onclick='gnuplot.toggle_plot("%s_plot_7")'>&#9318;</td>
    </tr>
    </table>
    </td>
  </tr>
  <tr style="border:0;">
    <td style="border:0;">
<canvas id="%s" width=600 height=400 tabindex="0">
    <div class='box'><h2>Your browser does not support the HTML 5 canvas element</h2></div>
</canvas>
    </td>
  </tr>
</table>
      <script>
          %s();
           $('body').on('contextmenu', '#%s', function(e){ return false; });
      </script>
      
            ''' % ((self._jsdir,)*4 + (imgdata,) + (self._jsdir,)*5 + (canvasname,)*12)

      #print(imgdata)
      img = showfunc(imgdata)
      if self.cleanup and not (self.suffix=='svg' and not uri):
         os.system("rm -f "+imgfile)
         # print(counter, end='\r')
         # print(counter, imgfile, os.path.exists(imgfile), os.system("lsof "+imgfile))
      else:
         print(imgfile)
      return img

# a default instance
gplot = Gplot()
ogplot = gplot.oplot


