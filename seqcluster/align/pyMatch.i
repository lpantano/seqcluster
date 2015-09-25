 %module pyMatch
 %{
 /* Put header files here or function declarations like below */
 extern int Match(char *strl, char *strs, int err, int err2);
 extern int MMatch(int len, char **ls, char *ss, int err, int err2);
 extern int Miraligner(char *lfn, char *sfn, char *rfn, int err, int err2);
 %}
 
 extern int Match(char *strl, char *strs, int err, int err2);
 extern int MMatch(int len, char **ls, char *ss, int err, int err2);
 extern int Miraligner(char *lfn, char *sfn, char *rfn, int err, int err2);
