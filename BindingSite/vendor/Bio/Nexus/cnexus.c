                                                       
  
                                                                       
                                                                             
                                                                             
           
  
           
  
                                                                         
                                                                         
                                                      
  
                   
   

#include <Python.h>
#include <string.h>

static PyObject * cnexus_scanfile(PyObject *self, PyObject *args)
{
    PyObject *cleaninput;
    const char *input;
    char *scanned, *scanned_start;
    char t, quotelevel;
    int speciallevel, commlevel;

    quotelevel=0;
    speciallevel=0;
    commlevel=0;

    if (!PyArg_ParseTuple(args, "s", &input))
        return NULL;
    if (!(scanned=PyMem_RawMalloc(strlen(input)+1)))
        PyErr_NoMemory();
    scanned_start=scanned;
    for(t=*input;(t=*input);input++)
    {
                                   
        if (!(commlevel || speciallevel) && t==quotelevel)
            quotelevel=0;
                                     
        else if (!quotelevel && !(commlevel || speciallevel) && (t=='\'' || t=='"'))
            quotelevel=t;
                                            
        else if (!quotelevel  && t=='[')
        {
                                            
                                                                           
                                                                           
                                                    
                               
              
            if ((*(input+1)=='&') && !(commlevel || speciallevel))
                speciallevel=1;
            else                       
                commlevel++;
        }
        else if (!quotelevel && t==']')
        {
                                                
            if (speciallevel)
                speciallevel=0;
            else
            {
                commlevel--;
                if (commlevel<0)                         
                {
                    PyMem_RawFree(scanned_start);
                    return Py_BuildValue("s","]");
                }
                continue;
            }
        }
        if (!commlevel)
        {
                                                                        
                                                                       
            if (t==';' && !(quotelevel || speciallevel))
                                                                              
                               
                *(scanned++)=7;
            else
                *(scanned++)=t;
        }
                                                                                       
                                                        
           
    }

    if (commlevel>0)
    {
                                
        PyMem_RawFree(scanned_start);
        return Py_BuildValue("s","[");
    }
    else
    {
        *scanned=0;                    
        cleaninput= Py_BuildValue("s",scanned_start);
        PyMem_RawFree(scanned_start);
        return cleaninput;
    }
}

static PyMethodDef cNexusMethods[] =
{
    {"scanfile",cnexus_scanfile,METH_VARARGS,"Scan file and deal with comments and quotes."},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "cnexus",
        NULL,
        -1,
        cNexusMethods,
        NULL,
        NULL,
        NULL,
        NULL
};

PyObject *
PyInit_cnexus(void)
{
    return PyModule_Create(&moduledef);
}
