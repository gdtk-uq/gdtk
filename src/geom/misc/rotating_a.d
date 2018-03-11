// rotating_a.d
// Copied directly from the plotutils manual with minimal adaption to D.
// PJ 2018-03-11
//
// Build with command:
// $ dmd -of=rotating_a rotating_a.d libplot.d -L=-lplot
//
import std.stdio;
import core.stdc.stdio: stdin, stdout, stderr;
import std.string: toStringz;
import core.time: dur;
import core.thread: Thread;
import libplot;

int main()
{
  int handle, angle = 0;

  /* set Plotter parameters */        
  pl_parampl ("BITMAPSIZE", cast(void*)toStringz("300x300")); 
  pl_parampl ("BG_COLOR", cast(void*)toStringz("blue")); /* background color for window */
  pl_parampl ("USE_DOUBLE_BUFFERING", cast(void*)toStringz("yes"));

  /* create an X Plotter with the specified parameters */
  handle = pl_newpl ("X", stdin, stdout, stderr);
  pl_selectpl (handle);

  /* open X Plotter, initialize coordinates, pen, and font */
  pl_openpl ();
  pl_fspace (0.0, 0.0, 1.0, 1.0);  /* use normalized coordinates */
  pl_pencolorname ("white");
  pl_ffontsize (1.0);
  pl_fontname ("NewCenturySchlbk-Roman");

  pl_fmove (.50,.50);           /* move to center */
  while (1)                     /* loop endlessly */
    {
      pl_erase ();
      pl_textangle (angle++);      /* set new rotation angle */
      pl_alabel ('c', 'c', "A");   /* draw a centered `A' */
      Thread.sleep(dur!("msecs")(100));
    }
  pl_closepl();                 /* close Plotter */

  pl_selectpl (0);              /* select default Plotter */
  pl_deletepl (handle);         /* delete Plotter we used */
  return 0;
}

