#include <stdio.h>
#include <string.h>


int Match(char *strl, char *strs, int err) {

	int lenl,lens;
	int error;
	int i,j;
	int pos = -1;
	
	lenl = strlen (strl);
	lens = strlen (strs);
	
	error = -1;
	
	i = 0;
	while (error < err && (i + lens) <= lenl) {
		j = 0;
		
		while (error < err && j < lens) {
			if (strl[i+j] != strs[j]) error++;
			j++;
		}
		
		if (error < err) {
			pos = i;
			break;
		}
		
		error = -1;
		
		i++;
	}
	
	return pos;
}

int main(int argc, char *argv) {
	printf("\nResult: %d", Match("Hola Como Estas", "Como", 0));
	printf("\nResult: %d", Match("Hola Como Estas", "Camo", 1));
	printf("\nResult: %d", Match("Hola Como Estas", "Estas", 0));
	printf("\n");
	
	return 1;
}