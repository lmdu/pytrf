#define PY_SSIZE_T_CLEAN
#include "tandem.h"

PyObject *find_short_tandem_repeats(PyObject *self, PyObject *args, PyObject *kwargs) {
	const char *seq;
	Py_ssize_t size;

	//minimum repeat array
	int minreps[7] = {0, 0, 0, 0, 0, 0, 0};

	// current slide position
	Py_ssize_t k, i;

	//motif length
	int m, j;

	//jump motif length
	int jump;

	//may be not a candiate repeat
	int nonrep;

	// tandem repeat start location
	Py_ssize_t start;

	// tandem repeat length
	int replen;

	// mono repeat len
	int mono_len;

	//tandem repeat repeats
	int repeats;

	//motif sequence
	char motif[7];

	static char *keywords[] = {"seq", "minreps", NULL};
	if (!PyArg_ParseTupleAndKeywords(args, kwargs, "s#(iiiiii)", keywords, &seq, &size, &minreps[1],
		&minreps[2], &minreps[3], &minreps[4], &minreps[5], &minreps[6])) {
		return NULL;
	}

	for (i = 1; i <= 6; ++i) {
		minreps[i] = minreps[i] * i;
	}

	PyObject *ssrs = PyList_New(0);
	PyObject *tmp;
	mono_len = 0;

	for (i = 0; i < size; ++i) {
		if (seq[i] == 78) {
			continue;
		}

		for (j = 1; j <= 6; ++j) {
			//check minimal repeat
			/*jump = j * (minreps[j] - 1);
			nonrep = 0;
			for (m = 0; m < j; m++) {
				k = i+m;
				if (seq[k] != seq[k+jump]) {
					nonrep = 1;
				}
			}

			if (nonrep) {
				continue;
			}*/

			start = i;

			while (seq[i] == seq[i+j]) {
				++i;
			}

			replen = i + j - start;

			printf("%c, %d, %d\n", seq[start], j, replen);

			if (replen >= minreps[j]) {
				memcpy(motif, seq+start, j);
				motif[j] = '\0';

				repeats = replen/j;

				i += j;
 
				if (motif[0] != 78 && i <= size) {
					tmp = Py_BuildValue("(siiKKi)", motif, j, repeats, start+1, i, replen);
					PyList_Append(ssrs, tmp);
					Py_DECREF(tmp);
				}

				j = 0;
			} else {
				i = start;
			}
		}
	}
	return ssrs;
}

