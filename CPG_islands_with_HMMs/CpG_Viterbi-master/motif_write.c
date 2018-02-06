//c script to map motif hits to cpg island

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

int main()
{
	FILE *fp;
	FILE *fp2 = fopen("results2.txt", "a+"); //w+?
	if (!fp2) {
		perror("fopen results2.txt");
		return 1;
	}
	FILE *fpout = fopen("results_out.txt", "w+"); //w+?
	if (!fpout) {
		perror("fopen results_out.txt");
		return 1;
	}
		

	//
	char *final[790];
	char *result;
	char *tempstr5 = "grep -w 'CpG Island ";
	char *tempstr6 = "' motif_hits.txt | wc -l";
	char eee[15];
	char *grepString;
	grepString = malloc(1025);
	int k;
	for (k = 1; k < 781; k++){
	//fp = popen("grep -w 'CpG Island 1' motif_hits2.txt | wc -l", "r");
		sprintf(eee, "%d", k);
		grepString = strcpy(grepString, tempstr5);
		grepString = strcat(grepString, eee);
		grepString = strcat(grepString, tempstr6);
		fp = popen(grepString, "r");
		//printf("%d\n", result);
		char buf[1024];
		while (fgets(buf, 1024, fp)){
			result = buf;
		}
		fclose(fp);
	//printf("%s", result);
	final[k] = strdup(result);
	}
	//while not eof
	char *temp;
	char *tempStr;
	char *Str = "CpG Island ";
	temp = malloc(1025);
	tempStr = malloc(1025);
	char buf2[1024];
	char e[15];
	int count = 0;

	while (fgets(buf2, 1024, fp2)){
		count += 1;
		sprintf(e, "%d", count);
		//e = (char)count;
		tempStr = strcpy(tempStr, Str);
		tempStr = strcat(tempStr, e);
		tempStr = strcat(tempStr, ":");
		if (strstr(buf2, tempStr)) {
			//printf("%s", buf2);
			size_t ln = strlen(buf2) - 1;
			if (buf2[ln] == '\n') {
				buf2[ln] = '\0';
			}
			temp = strcpy(temp, buf2);
			temp = strcat(temp, "; ");
			temp = strcat(temp, final[count]);
			//temp = strcat(temp, "\n");
			fprintf(fpout, "%s", temp);
		}
	}

	fclose(fp);
	fclose(fp2);

	return 0;
}
