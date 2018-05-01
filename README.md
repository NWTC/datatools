# datatools

Repository hosting python code developed at NWTC for data processing, 
for example, input/output of commonly used file formats.

Code-dependent outputs are organized into submodules:

* datatools.SOWFA
* datatools.FAST
* etc.

Helper scripts to be run from the command line that utilize these tools
will live in datatools/utilities

Acknowledgment
==============
If you use this code as part of any published research, please
acknowledge the ``NWTC Open Source Wind Research Repositories
(OSWRR)``. 

Guidelines
==============
        All added code must be documented with python docstrings as shown in the example below.
        Ideally, you should also provide a minimal working example to make your code more easily usable by others.
                

        def function_name(arg1, arg2):
        """
        Summary line.

        Extended description of function.

        Parameters
        ----------
        arg1 : int
                Description of arg1
        arg2 : str
                Description of arg2

        Returns
        -------
        int
                Description of return value


        @author: your_github_username
        """
        return some_variable

