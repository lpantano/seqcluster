 %module pyMatch
 %{
 /* Put header files here or function declarations like below */
 extern int Match(char *strl, char *strs, int err);
 %}
 
 extern int Match(char *strl, char *strs, int err);
