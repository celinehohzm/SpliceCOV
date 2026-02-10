#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// Global parameters
static double lowcov = 10;
static double splicenoise = 0.01;
static double percnoise = 0.2;
static double highnoise = 0.005;
static double smallcov = 50;
static double highcov = 100;
static int delta_param = 5;
static int win = 150;
static int smallwin = 25;

// Data structures
typedef struct {
    int start;
    int end;
    double cov;
} CovgEntry;

typedef struct {
    char chrname[256];
    int start;
    int end;
    double cov;
    char strand;
    double ps;
} JuncEntry;

typedef struct {
    int pos;
    double perc;
    double covdiff;
} DropEntry;

typedef struct {
    char type[8];  // "tstart", "tend", "jstart", "jend"
    int pos;
    int *indices;
    int num_indices;
    double change_perc;
    double pos_cov;
    double cov_to_next;
} RecordEntry;

typedef struct {
    int pos;
    double perc;
    double cov;
    int active;
} TmpEntry;

// Dynamic arrays
static CovgEntry *covg = NULL;
static int covg_size = 0;
static int covg_cap = 0;

static JuncEntry *junc = NULL;
static int junc_size = 0;
static int junc_cap = 0;

static JuncEntry *unprocjunc = NULL;
static int unprocjunc_size = 0;
static int unprocjunc_cap = 0;

static DropEntry *drop_arr = NULL;
static int drop_size = 0;
static int drop_cap = 0;

static RecordEntry *record = NULL;
static int record_size = 0;
static int record_cap = 0;

static int *jend = NULL;
static int jend_size = 0;

// Function prototypes
static void push_covg(int start, int end, double cov);
static void push_junc(JuncEntry *arr, int *size, int *cap, const char *chrname, int start, int end, double cov, char strand, double ps);
static void push_drop(int pos, double perc, double covdiff);
static void push_record(const char *type, int pos, int *indices, int num_indices, double change_perc, double pos_cov, double cov_to_next);
static void clear_covg(void);
static void clear_junc(void);
static void clear_drop(void);
static void clear_record(void);
static void clear_unprocjunc(void);

static int process_bundle(const char *chr, int bundleno);
static int process_records(const char *chr, int bundleno);
static int print_small_bundle(const char *chr, int bundleno, int s, int e);
static double get_next_cov(int i, int e);
static void get_record(int *ib, int *id, int *js, int *je, int *prevpos, int nb, int nd, int nj, const char *chr);
static int less_than(int n1, int n2);
static void get_drop(int si, int se, int nb, int *js_out, int *je_out, int js, int je, int nj);
static void compute_perc(int l, int r, int i, double *cov, double *adjs, double *adje, double *percl, double *percr, double *avgl, double *avgr);
static double get_cov(int start, int end, int *si, int nb);
static int add_procjunc_to_bundle(int *bundleend, const char *chr);
static int add_junc_to_bundle(const char *chr, int bundleend, FILE *fJ);
static void process_junctions(int nj);
static int equal_strand(char s1, char s2);
static void sort_jend(int nj);
static int compare_jend(const void *a, const void *b);

// Helper functions for dynamic arrays
static void push_covg(int start, int end, double cov) {
    if (covg_size >= covg_cap) {
        covg_cap = covg_cap ? covg_cap * 2 : 1024;
        covg = realloc(covg, covg_cap * sizeof(CovgEntry));
    }
    covg[covg_size].start = start;
    covg[covg_size].end = end;
    covg[covg_size].cov = cov;
    covg_size++;
}

static void push_junc(JuncEntry *arr, int *size, int *cap, const char *chrname, int start, int end, double cov, char strand, double ps) {
    if (*size >= *cap) {
        *cap = *cap ? *cap * 2 : 1024;
        if (arr == junc) {
            junc = realloc(junc, *cap * sizeof(JuncEntry));
            arr = junc;
        } else {
            unprocjunc = realloc(unprocjunc, *cap * sizeof(JuncEntry));
            arr = unprocjunc;
        }
    }
    strncpy(arr[*size].chrname, chrname, 255);
    arr[*size].chrname[255] = '\0';
    arr[*size].start = start;
    arr[*size].end = end;
    arr[*size].cov = cov;
    arr[*size].strand = strand;
    arr[*size].ps = ps;
    (*size)++;
}

static void push_junc_main(const char *chrname, int start, int end, double cov, char strand, double ps) {
    if (junc_size >= junc_cap) {
        junc_cap = junc_cap ? junc_cap * 2 : 1024;
        junc = realloc(junc, junc_cap * sizeof(JuncEntry));
    }
    strncpy(junc[junc_size].chrname, chrname, 255);
    junc[junc_size].chrname[255] = '\0';
    junc[junc_size].start = start;
    junc[junc_size].end = end;
    junc[junc_size].cov = cov;
    junc[junc_size].strand = strand;
    junc[junc_size].ps = ps;
    junc_size++;
}

static void push_unprocjunc(const char *chrname, int start, int end, double cov, char strand, double ps) {
    if (unprocjunc_size >= unprocjunc_cap) {
        unprocjunc_cap = unprocjunc_cap ? unprocjunc_cap * 2 : 1024;
        unprocjunc = realloc(unprocjunc, unprocjunc_cap * sizeof(JuncEntry));
    }
    strncpy(unprocjunc[unprocjunc_size].chrname, chrname, 255);
    unprocjunc[unprocjunc_size].chrname[255] = '\0';
    unprocjunc[unprocjunc_size].start = start;
    unprocjunc[unprocjunc_size].end = end;
    unprocjunc[unprocjunc_size].cov = cov;
    unprocjunc[unprocjunc_size].strand = strand;
    unprocjunc[unprocjunc_size].ps = ps;
    unprocjunc_size++;
}

static void push_drop(int pos, double perc, double covdiff) {
    if (drop_size >= drop_cap) {
        drop_cap = drop_cap ? drop_cap * 2 : 1024;
        drop_arr = realloc(drop_arr, drop_cap * sizeof(DropEntry));
    }
    drop_arr[drop_size].pos = pos;
    drop_arr[drop_size].perc = perc;
    drop_arr[drop_size].covdiff = covdiff;
    drop_size++;
}

static void push_record(const char *type, int pos, int *indices, int num_indices, double change_perc, double pos_cov, double cov_to_next) {
    if (record_size >= record_cap) {
        record_cap = record_cap ? record_cap * 2 : 1024;
        record = realloc(record, record_cap * sizeof(RecordEntry));
    }
    strncpy(record[record_size].type, type, 7);
    record[record_size].type[7] = '\0';
    record[record_size].pos = pos;
    record[record_size].indices = malloc(num_indices * sizeof(int));
    memcpy(record[record_size].indices, indices, num_indices * sizeof(int));
    record[record_size].num_indices = num_indices;
    record[record_size].change_perc = change_perc;
    record[record_size].pos_cov = pos_cov;
    record[record_size].cov_to_next = cov_to_next;
    record_size++;
}

static void clear_covg(void) {
    covg_size = 0;
}

static void clear_junc(void) {
    junc_size = 0;
}

static void clear_drop(void) {
    drop_size = 0;
}

static void clear_record(void) {
    for (int i = 0; i < record_size; i++) {
        free(record[record_size].indices);
    }
    record_size = 0;
}

static void clear_unprocjunc(void) {
    unprocjunc_size = 0;
}

// Comparison function for sorting jend
static int compare_jend(const void *a, const void *b) {
    int ia = *(const int *)a;
    int ib = *(const int *)b;
    return junc[ia].end - junc[ib].end;
}

static void sort_jend(int nj) {
    if (jend) free(jend);
    jend = malloc(nj * sizeof(int));
    jend_size = nj;
    for (int i = 0; i < nj; i++) {
        jend[i] = i;
    }
    qsort(jend, nj, sizeof(int), compare_jend);
}

static int equal_strand(char s1, char s2) {
    if (s1 == s2) return 1;
    if (s1 == '.' || s2 == '.') return 1;
    return 0;
}

static int less_than(int n1, int n2) {
    if (n1 == 0) return 0;
    if (n2 == 0) return 1;
    if (n1 < n2) return 1;
    return 0;
}

static double get_cov(int start, int end, int *si, int nb) {
    double cov_sum = 0;
    
    while (*si < nb && start > covg[*si].end) {
        (*si)++;
    }
    
    if (*si == nb) return cov_sum;
    
    if (start < covg[*si].start) {
        start = covg[*si].start;
    }
    
    while (end > covg[*si].end) {
        cov_sum += (covg[*si].end - start + 1) * covg[*si].cov;
        (*si)++;
        if (*si == nb) return cov_sum;
        start = covg[*si].start;
    }
    
    if (end < start) return cov_sum;
    
    cov_sum += (end - start + 1) * covg[*si].cov;
    
    return cov_sum;
}

static void process_junctions(int nj) {
    typedef struct {
        char key[512];
        int idx;
    } ActiveEntry;
    
    ActiveEntry *active = malloc(nj * sizeof(ActiveEntry));
    int active_size = 0;
    
    int *mark = calloc(nj, sizeof(int));
    
    int js = 0;
    int je = 0;
    
    while (js < nj) {
        // Remove junctions that ended before current start
        while (je < nj && junc[jend[je]].end < junc[js].start) {
            char key[512];
            snprintf(key, sizeof(key), "%s:%d:%c", junc[jend[je]].chrname, junc[jend[je]].start, junc[jend[je]].strand);
            
            for (int i = 0; i < active_size; i++) {
                if (strcmp(active[i].key, key) == 0) {
                    // Remove by shifting
                    for (int j = i; j < active_size - 1; j++) {
                        active[j] = active[j + 1];
                    }
                    active_size--;
                    break;
                }
            }
            je++;
        }
        
        // Check against active junctions
        for (int i = 0; i < active_size; i++) {
            int ai = active[i].idx;
            if (equal_strand(junc[ai].strand, junc[js].strand)) {
                if (junc[js].cov > junc[ai].cov) {
                    if (junc[js].cov < smallcov) {
                        mark[ai] = 1;
                    }
                } else if (junc[js].cov < junc[ai].cov) {
                    if (junc[ai].cov < smallcov) {
                        mark[js] = 1;
                    }
                }
            }
        }
        
        // Add current junction to active
        char key[512];
        snprintf(key, sizeof(key), "%s:%d:%c", junc[js].chrname, junc[js].start, junc[js].strand);
        strcpy(active[active_size].key, key);
        active[active_size].idx = js;
        active_size++;
        
        js++;
    }
    
    // Apply marks
    for (int j = 0; j < nj; j++) {
        if (mark[j]) {
            junc[j].cov = 0;
        }
    }
    
    free(active);
    free(mark);
}

static void compute_perc(int l, int r, int i, double *cov, double *adjs, double *adje, 
                         double *percl, double *percr, double *avgl, double *avgr) {
    *percl = 1;
    *percr = 1;
    *avgl = 0;
    *avgr = 0;
    
    double sumleft = cov[i - 1] - cov[l - 1];
    double sumlefta = adjs[i - 1] - adjs[l - 1];
    int k = (i - 2) % win;
    int kw = i - 1 - k;
    if (kw < i && kw >= l) {
        sumleft += cov[kw - 1];
        sumlefta += adjs[kw - 1];
    }
    sumlefta = sumleft - sumlefta;
    
    double sumright = cov[r] - cov[i - 1];
    double sumrighta = adje[r] - adje[i - 1];
    k = (r - 1) % win;
    kw = r - k;
    if (kw - 1 < r && kw >= i) {
        sumright += cov[kw - 1];
        sumrighta += adje[kw - 1];
    }
    sumrighta = sumright - sumrighta;
    
    if (sumlefta > sumright) {
        *avgl = (sumlefta - sumright) / (i - l);
        *percl = sumright / sumlefta;
    }
    
    if (sumrighta > sumleft) {
        *avgr = (sumrighta - sumleft) / (i - l);
        *percr = sumleft / sumrighta;
    }
}

static void get_drop(int si, int se, int nb, int *js_out, int *je_out, int js, int je, int nj) {
    int start = covg[si].start;
    int end = covg[se].end;
    
    while (js < nj && junc[js].start < start) js++;
    while (je < nj && junc[jend[je]].end < start) je++;
    
    push_drop(start, 0, covg[si].cov);
    
    int len = end - start + 1;
    if (len >= win + smallwin) {
        // Allocate arrays
        double *c = calloc(len + win + 10, sizeof(double));
        double *as = calloc(len + win + 10, sizeof(double));
        double *ae = calloc(len + win + 10, sizeof(double));
        
        int *jp = NULL;
        int jp_size = 0;
        int jp_cap = 0;
        
        c[0] = 0;
        as[0] = 0;
        ae[0] = 0;
        
        int curr_si = si;
        
        for (int i = 0; i < len; i++) {
            int istart = i + start;
            while (istart > covg[curr_si].end) curr_si++;
            int i1 = i + 1;
            c[i1] = covg[curr_si].cov;
            int r_mod = i % win;
            if (r_mod) {
                c[i1] += c[i];
                as[i1] += as[i];
                ae[i1] += ae[i];
            }
            
            int isjunc = 0;
            if (js < nj && junc[js].start == istart) {
                double jcov = junc[js].cov;
                js++;
                while (js < nj && junc[js].start == istart) {
                    jcov += junc[js].cov;
                    js++;
                }
                
                for (int j = 1; j < r_mod + 2; j++) {
                    as[i - r_mod + j] += jcov * j;
                }
                
                if (i > win) {
                    for (int j = 1; j < win - r_mod; j++) {
                        as[i - win + j + 1] += jcov * j;
                    }
                }
                
                // Push to jp
                if (jp_size >= jp_cap) {
                    jp_cap = jp_cap ? jp_cap * 2 : 256;
                    jp = realloc(jp, jp_cap * sizeof(int));
                }
                jp[jp_size++] = i1;
                isjunc = 1;
            }
            
            if (je < nj && junc[jend[je]].end == istart) {
                double jcov = junc[jend[je]].cov;
                je++;
                while (je < nj && junc[jend[je]].end == istart) {
                    jcov += junc[jend[je]].cov;
                    je++;
                }
                
                for (int j = 0; j < win; j++) {
                    ae[i1 + j] += jcov;
                }
                if (!isjunc) {
                    if (jp_size >= jp_cap) {
                        jp_cap = jp_cap ? jp_cap * 2 : 256;
                        jp = realloc(jp, jp_cap * sizeof(int));
                    }
                    jp[jp_size++] = i1;
                }
            }
        }
        
        int j_idx = 0;
        int i = smallwin + 1;
        
        TmpEntry *tmps = NULL;
        int tmps_size = 0;
        int tmps_cap = 0;
        
        TmpEntry *tmpe = NULL;
        int tmpe_size = 0;
        int tmpe_cap = 0;
        
        int maxs = 0;
        int maxe = 0;
        
        while (i < len - smallwin) {
            int l = i - win;
            if (l < 1) l = 1;
            
            int r_val = 2 * i - l - 1;
            if (r_val > len - 1) {
                r_val = len - 1;
                l = 2 * i - r_val - 1;
            }
            
            double minpercl, minpercr, minavgl, minavgr;
            compute_perc(l, r_val, i, c, as, ae, &minpercl, &minpercr, &minavgl, &minavgr);
            
            while (j_idx < jp_size && jp[j_idx] <= l) j_idx++;
            
            int k = j_idx;
            while (k < jp_size && jp[k] < r_val) {
                if (i - jp[k] > smallwin) {
                    double percl, percr, avgl, avgr;
                    compute_perc(jp[k], 2 * i - jp[k] - 1, i, c, as, ae, &percl, &percr, &avgl, &avgr);
                    if (percl < minpercl) { minpercl = percl; minavgl = avgl; }
                    if (percr < minpercr) { minpercr = percr; minavgr = avgr; }
                } else if (jp[k] - i + 1 > smallwin) {
                    double percl, percr, avgl, avgr;
                    compute_perc(2 * i - jp[k] - 1, jp[k], i, c, as, ae, &percl, &percr, &avgl, &avgr);
                    if (percl < minpercl) { minpercl = percl; minavgl = avgl; }
                    if (percr < minpercr) { minpercr = percr; minavgr = avgr; }
                }
                k++;
            }
            
            int plus = 0;
            if (minpercr < percnoise) {
                if (tmps_size > 0) {
                    if (tmps_size >= tmps_cap) {
                        tmps_cap = tmps_cap ? tmps_cap * 2 : 256;
                        tmps = realloc(tmps, tmps_cap * sizeof(TmpEntry));
                    }
                    tmps[tmps_size].pos = i;
                    tmps[tmps_size].perc = minpercr;
                    tmps[tmps_size].cov = minavgr;
                    tmps[tmps_size].active = 1;
                    tmps_size++;
                    
                    if (i - tmps[maxs].pos > win) {
                        if (tmps[tmps_size - 2].perc >= minpercr || i - tmps[tmps_size - 2].pos > win) {
                            maxs = tmps_size - 1;
                            k = maxs - 1;
                            while (k >= 0 && i - tmps[k].pos <= win) {
                                if (minpercr > tmps[k].perc) {
                                    maxs = k;
                                    tmps[tmps_size - 1].active = 0;
                                } else {
                                    tmps[k].active = 0;
                                }
                                k--;
                            }
                        } else {
                            tmps[tmps_size - 1].active = 0;
                        }
                    } else if (minpercr < tmps[maxs].perc) {
                        tmps[maxs].active = 0;
                        maxs = tmps_size - 1;
                    } else {
                        tmps[tmps_size - 1].active = 0;
                    }
                } else {
                    if (tmps_size >= tmps_cap) {
                        tmps_cap = tmps_cap ? tmps_cap * 2 : 256;
                        tmps = realloc(tmps, tmps_cap * sizeof(TmpEntry));
                    }
                    tmps[tmps_size].pos = i;
                    tmps[tmps_size].perc = minpercr;
                    tmps[tmps_size].cov = minavgr;
                    tmps[tmps_size].active = 1;
                    tmps_size++;
                }
                plus = 1;
            }
            
            if (minpercl < percnoise) {
                if (tmpe_size > 0) {
                    if (tmpe_size >= tmpe_cap) {
                        tmpe_cap = tmpe_cap ? tmpe_cap * 2 : 256;
                        tmpe = realloc(tmpe, tmpe_cap * sizeof(TmpEntry));
                    }
                    tmpe[tmpe_size].pos = i;
                    tmpe[tmpe_size].perc = minpercl;
                    tmpe[tmpe_size].cov = minavgl;
                    tmpe[tmpe_size].active = 1;
                    tmpe_size++;
                    
                    if (i - tmpe[maxe].pos > win) {
                        if (tmpe[tmpe_size - 2].perc >= minpercl || i - tmpe[tmpe_size - 2].pos > win) {
                            maxe = tmpe_size - 1;
                            k = maxe - 1;
                            while (k >= 0 && i - tmpe[k].pos <= win) {
                                if (minpercl > tmpe[k].perc) {
                                    maxe = k;
                                    tmpe[tmpe_size - 1].active = 0;
                                } else {
                                    tmpe[k].active = 0;
                                }
                                k--;
                            }
                        } else {
                            tmpe[tmpe_size - 1].active = 0;
                        }
                    } else if (minpercl < tmpe[maxe].perc) {
                        tmpe[maxe].active = 0;
                        maxe = tmpe_size - 1;
                    } else {
                        tmpe[tmpe_size - 1].active = 0;
                    }
                } else {
                    if (tmpe_size >= tmpe_cap) {
                        tmpe_cap = tmpe_cap ? tmpe_cap * 2 : 256;
                        tmpe = realloc(tmpe, tmpe_cap * sizeof(TmpEntry));
                    }
                    tmpe[tmpe_size].pos = i;
                    tmpe[tmpe_size].perc = minpercl;
                    tmpe[tmpe_size].cov = minavgl;
                    tmpe[tmpe_size].active = 1;
                    tmpe_size++;
                }
                plus = 1;
            }
            
            if (!plus && minpercl > 0.5 && minpercr > 0.5) {
                plus = delta_param;
            } else if (!plus) {
                plus = 1;
            }
            
            i += plus;
        }
        
        int ns = tmps_size;
        int ne = tmpe_size;
        
        int s_idx = 0;
        int e_idx = 0;
        
        while (s_idx < ns && e_idx < ne) {
            if (tmpe[e_idx].pos < tmps[s_idx].pos) {
                if (tmpe[e_idx].active) {
                    push_drop(-(tmpe[e_idx].pos + start - 2), tmpe[e_idx].perc, tmpe[e_idx].cov);
                }
                e_idx++;
            } else {
                if (tmps[s_idx].active) {
                    if (tmpe[e_idx].active && tmpe[e_idx].pos == tmps[s_idx].pos) {
                        exit(1);
                    }
                    push_drop(tmps[s_idx].pos + start - 1, tmps[s_idx].perc, tmps[s_idx].cov);
                }
                s_idx++;
            }
        }
        while (s_idx < ns) {
            if (tmps[s_idx].active) {
                push_drop(tmps[s_idx].pos + start - 1, tmps[s_idx].perc, tmps[s_idx].cov);
            }
            s_idx++;
        }
        while (e_idx < ne) {
            if (tmpe[e_idx].active) {
                push_drop(-(tmpe[e_idx].pos + start - 2), tmpe[e_idx].perc, tmpe[e_idx].cov);
            }
            e_idx++;
        }
        
        free(c);
        free(as);
        free(ae);
        free(jp);
        free(tmps);
        free(tmpe);
    }
    
    if (js > 0) js--;
    while (js < nj && junc[js].start < end) js++;
    
    push_drop(-end, 0, covg[se].cov);
    
    *js_out = js;
    *je_out = je;
}

static void get_record(int *ib, int *id, int *js, int *je, int *prevpos, int nb, int nd, int nj, const char *chr) {
    int nextd = 0;
    int nextjs = 0;
    int nextje = 0;
    
    if (*id < nd) {
        nextd = abs(drop_arr[*id].pos);
    }
    
    while (*js < nj) {
        if (junc[*js].cov > 0) {
            nextjs = junc[*js].start;
            break;
        }
        (*js)++;
    }
    
    while (*je < nj) {
        if (junc[jend[*je]].cov > 0) {
            nextje = junc[jend[*je]].end;
            break;
        }
        (*je)++;
    }
    
    if (!nextd && !nextjs && !nextje) {
        if (*id == nd && *js == nj && *je == nj) return;
        exit(1);
    }
    
    if (less_than(nextd, nextjs)) {
        if (less_than(nextd, nextje)) {
            // start/stop is smallest
            if (nextd > *prevpos + 1) {
                int tmpib = *ib;
                double avgcov = get_cov(*prevpos + 1, nextd - 1, &tmpib, nb);
                record[record_size - 1].cov_to_next += avgcov;
            }
            
            int tmpib = *ib;
            double pos_cov = get_cov(nextd, nextd, &tmpib, nb);
            *ib = tmpib;
            
            while (*id < nd && abs(drop_arr[*id].pos) == nextd) {
                const char *type = "tstart";
                if (drop_arr[*id].pos < 0) type = "tend";
                
                // Check if should delete previous junctions
                if (drop_arr[*id].perc == 0 && *prevpos > nextd - delta_param) {
                    if (strcmp(type, "tstart") == 0) {
                        if (record[record_size - 1].pos == nextd && strcmp(record[record_size - 1].type, "jstart") == 0) {
                            for (int j = 0; j < record[record_size - 1].num_indices; j++) {
                                junc[record[record_size - 1].indices[j]].cov = 0;
                            }
                        }
                    } else {
                        int nr = record_size;
                        int i = nr - 1;
                        while (i >= 0 && record[i].pos > nextd - delta_param) {
                            if (strcmp(record[i].type, "jend") == 0) {
                                for (int j = 0; j < record[i].num_indices; j++) {
                                    junc[record[i].indices[j]].cov = 0;
                                }
                            }
                            i--;
                        }
                    }
                }
                
                int idx = *id;
                push_record(type, nextd, &idx, 1, drop_arr[*id].perc, pos_cov, 0);
                
                (*id)++;
            }
            
            *prevpos = nextd;
        } else if (nextje) {
            // nextje is smallest
            int *tmpr = NULL;
            int tmpr_size = 0;
            int tmpr_cap = 0;
            double count = 0;
            int present = 0;
            
            while (*je < nj && junc[jend[*je]].end == nextje) {
                if (junc[jend[*je]].cov > 0) {
                    if (tmpr_size >= tmpr_cap) {
                        tmpr_cap = tmpr_cap ? tmpr_cap * 2 : 16;
                        tmpr = realloc(tmpr, tmpr_cap * sizeof(int));
                    }
                    tmpr[tmpr_size++] = jend[*je];
                    count += junc[jend[*je]].cov;
                }
                (*je)++;
            }
            
            if (count > 0) {
                present = 1;
                
                int tmpib = *ib;
                int leftstart = nextje - delta_param;
                if (leftstart < covg[0].start) {
                    leftstart = covg[0].start;
                    tmpib = 0;
                } else {
                    while (covg[tmpib].start > leftstart) tmpib--;
                }
                double leftcov = get_cov(leftstart, nextje - 1, &tmpib, nb);
                
                int rightend = nextje + delta_param - 1;
                if (rightend > covg[nb - 1].end) rightend = covg[nb - 1].end;
                double rightcov = get_cov(nextje, rightend, &tmpib, nb);
                
                if (leftcov < rightcov) {
                    double prevcount = 0;
                    int prevjend = 0;
                    if (*prevpos > nextje - 2) {
                        int nr = record_size;
                        int i = nr - 1;
                        while (i >= 0 && record[i].pos > nextje - 2) {
                            if (record[i].pos == nextje - 1 && strcmp(record[i].type, "jend") == 0) {
                                prevjend = i;
                                for (int j = 0; j < record[i].num_indices; j++) {
                                    prevcount += junc[record[i].indices[j]].cov;
                                }
                            }
                            i--;
                        }
                    }
                    
                    if (prevcount < count) {
                        if (prevcount > 0) {
                            for (int j = 0; j < record[prevjend].num_indices; j++) {
                                junc[record[prevjend].indices[j]].cov = 0;
                            }
                        }
                        
                        if (nextje > *prevpos + 1) {
                            int tmp_ib = *ib;
                            double avgcov = get_cov(*prevpos + 1, nextje - 1, &tmp_ib, nb);
                            record[record_size - 1].cov_to_next += avgcov;
                        }
                        
                        tmpib = *ib;
                        double pos_cov = get_cov(nextje, nextje, &tmpib, nb);
                        *ib = tmpib;
                        
                        push_record("jend", nextje, tmpr, tmpr_size, leftcov / rightcov, pos_cov, 0);
                        
                        *prevpos = nextje;
                    } else {
                        count = 0;
                    }
                } else {
                    count = 0;
                }
            }
            
            if (present && count == 0) {
                for (int j = 0; j < tmpr_size; j++) {
                    junc[tmpr[j]].cov = 0;
                }
            }
            
            free(tmpr);
        }
    } else {
        if (less_than(nextjs, nextje)) {
            // nextjs is smallest
            int *tmpr = NULL;
            int tmpr_size = 0;
            int tmpr_cap = 0;
            double count = 0;
            int present = 0;
            
            while (*js < nj && junc[*js].start == nextjs) {
                if (junc[*js].cov > 0) {
                    if (tmpr_size >= tmpr_cap) {
                        tmpr_cap = tmpr_cap ? tmpr_cap * 2 : 16;
                        tmpr = realloc(tmpr, tmpr_cap * sizeof(int));
                    }
                    tmpr[tmpr_size++] = *js;
                    count += junc[*js].cov;
                }
                (*js)++;
            }
            
            if (count > 0) {
                present = 1;
                
                double prevcount = 0;
                int prevjstart = 0;
                if (*prevpos > nextjs - delta_param) {
                    int nr = record_size;
                    int i = nr - 1;
                    while (i >= 0 && record[i].pos > nextjs - delta_param) {
                        if ((record[i].change_perc == 0 && strcmp(record[i].type, "tstart") == 0) ||
                            (record[i].pos == nextjs && strcmp(record[i].type, "jend") == 0)) {
                            count = 0;
                            break;
                        } else if (record[i].pos == nextjs - 1 && strcmp(record[i].type, "jstart") == 0) {
                            prevjstart = i;
                            for (int j = 0; j < record[i].num_indices; j++) {
                                prevcount += junc[record[i].indices[j]].cov;
                            }
                        }
                        i--;
                    }
                }
                
                if (prevcount < count && count > 0) {
                    int tmpib = *ib;
                    int leftstart = nextjs - delta_param + 1;
                    if (leftstart < covg[0].start) {
                        leftstart = covg[0].start;
                        tmpib = 0;
                    } else {
                        while (covg[tmpib].start > leftstart) tmpib--;
                    }
                    double leftcov = get_cov(leftstart, nextjs, &tmpib, nb);
                    
                    int rightend = nextjs + delta_param;
                    if (rightend > covg[nb - 1].end) rightend = covg[nb - 1].end;
                    double rightcov = get_cov(nextjs + 1, rightend, &tmpib, nb);
                    
                    if (leftcov > rightcov) {
                        if (prevcount > 0) {
                            for (int j = 0; j < record[prevjstart].num_indices; j++) {
                                junc[record[prevjstart].indices[j]].cov = 0;
                            }
                        }
                        
                        if (nextjs > *prevpos + 1) {
                            int tmp_ib = *ib;
                            double avgcov = get_cov(*prevpos + 1, nextjs - 1, &tmp_ib, nb);
                            record[record_size - 1].cov_to_next += avgcov;
                        }
                        
                        tmpib = *ib;
                        double pos_cov = get_cov(nextjs, nextjs, &tmpib, nb);
                        *ib = tmpib;
                        
                        push_record("jstart", nextjs, tmpr, tmpr_size, rightcov / leftcov, pos_cov, 0);
                        
                        *prevpos = nextjs;
                    } else {
                        count = 0;
                    }
                } else {
                    count = 0;
                }
            }
            
            if (present && count == 0) {
                for (int j = 0; j < tmpr_size; j++) {
                    junc[tmpr[j]].cov = 0;
                }
            }
            
            free(tmpr);
        } else if (nextje) {
            // nextje is smallest (same as above case)
            int *tmpr = NULL;
            int tmpr_size = 0;
            int tmpr_cap = 0;
            double count = 0;
            int present = 0;
            
            while (*je < nj && junc[jend[*je]].end == nextje) {
                if (junc[jend[*je]].cov > 0) {
                    if (tmpr_size >= tmpr_cap) {
                        tmpr_cap = tmpr_cap ? tmpr_cap * 2 : 16;
                        tmpr = realloc(tmpr, tmpr_cap * sizeof(int));
                    }
                    tmpr[tmpr_size++] = jend[*je];
                    count += junc[jend[*je]].cov;
                }
                (*je)++;
            }
            
            if (count > 0) {
                present = 1;
                
                int tmpib = *ib;
                int leftstart = nextje - delta_param;
                if (leftstart < covg[0].start) {
                    leftstart = covg[0].start;
                    tmpib = 0;
                } else {
                    while (covg[tmpib].start > leftstart) tmpib--;
                }
                double leftcov = get_cov(leftstart, nextje - 1, &tmpib, nb);
                
                int rightend = nextje + delta_param - 1;
                if (rightend > covg[nb - 1].end) rightend = covg[nb - 1].end;
                double rightcov = get_cov(nextje, rightend, &tmpib, nb);
                
                if (leftcov < rightcov) {
                    double prevcount = 0;
                    int prevjend = 0;
                    if (*prevpos > nextje - 2) {
                        int nr = record_size;
                        int i = nr - 1;
                        while (i >= 0 && record[i].pos > nextje - 2) {
                            if (record[i].pos == nextje - 1 && strcmp(record[i].type, "jend") == 0) {
                                prevjend = i;
                                for (int j = 0; j < record[i].num_indices; j++) {
                                    prevcount += junc[record[i].indices[j]].cov;
                                }
                            }
                            i--;
                        }
                    }
                    
                    if (prevcount < count) {
                        if (prevcount > 0) {
                            for (int j = 0; j < record[prevjend].num_indices; j++) {
                                junc[record[prevjend].indices[j]].cov = 0;
                            }
                        }
                        
                        if (nextje > *prevpos + 1) {
                            int tmp_ib = *ib;
                            double avgcov = get_cov(*prevpos + 1, nextje - 1, &tmp_ib, nb);
                            record[record_size - 1].cov_to_next += avgcov;
                        }
                        
                        tmpib = *ib;
                        double pos_cov = get_cov(nextje, nextje, &tmpib, nb);
                        *ib = tmpib;
                        
                        push_record("jend", nextje, tmpr, tmpr_size, leftcov / rightcov, pos_cov, 0);
                        
                        *prevpos = nextje;
                    } else {
                        count = 0;
                    }
                } else {
                    count = 0;
                }
            }
            
            if (present && count == 0) {
                for (int j = 0; j < tmpr_size; j++) {
                    junc[tmpr[j]].cov = 0;
                }
            }
            
            free(tmpr);
        }
    }
}

static double get_next_cov(int i, int e) {
    double cov_val = record[i].cov_to_next;
    int start = record[i].pos;
    
    if (strcmp(record[i].type, "tstart") == 0 || strcmp(record[i].type, "jend") == 0) {
        cov_val += record[i].pos_cov;
    } else {
        start++;
    }
    
    i++;
    while (i < e && record[i].pos == 0) {
        cov_val += record[i].pos_cov + record[i].cov_to_next;
        i++;
    }
    
    int len_val = 0;
    
    if (i < e) {
        if (strcmp(record[i].type, "tend") == 0 || strcmp(record[i].type, "jstart") == 0) {
            cov_val += record[i].pos_cov;
            len_val += 1;
        }
    } else {
        return 0;
    }
    
    len_val += record[i].pos - start;
    
    if (len_val > 0) {
        cov_val /= len_val;
    }
    
    return cov_val;
}

static int print_small_bundle(const char *chr, int bundleno, int s, int e) {
    int b = s;
    double sum = 0;
    double sumb = 0;
    int found = 0;
    int nl = 0;
    int reale = e - 1;
    
    for (int i = s; i < e; i++) {
        if (record[i].pos == 0) {
            sum += record[i].pos_cov + record[i].cov_to_next;
            if (b >= 0) {
                sumb += record[i].pos_cov + record[i].cov_to_next;
            }
        } else if (strcmp(record[i].type, "tend") == 0) {
            sum += record[i].pos_cov + record[i].cov_to_next;
            sumb += record[i].pos_cov + record[i].cov_to_next;
            nl++;
            if (record[i].change_perc == 0) {
                int len_val = record[i].pos - record[b].pos + 1;
                if (!found && (len_val < win || sumb / len_val < lowcov)) {
                    for (int j = b; j <= i; j++) {
                        if (record[j].pos) {
                            record[j].pos = 0;
                            nl--;
                        }
                    }
                } else {
                    reale = i;
                }
                b = -1;
                sumb = 0;
                found = 0;
            }
        } else if (strcmp(record[i].type, "tstart") == 0) {
            sum += record[i].pos_cov + record[i].cov_to_next;
            sumb += record[i].pos_cov + record[i].cov_to_next;
            nl++;
            if (record[i].change_perc == 0) {
                b = i;
                found = 0;
            }
        } else {
            sum += record[i].pos_cov + record[i].cov_to_next;
            sumb += record[i].pos_cov + record[i].cov_to_next;
            nl++;
            found = 1;
        }
    }
    
    if (nl > 0) {
        double avg = sum / (record[reale].pos - record[s].pos + 1);
        
        printf("bundle\t%s\t%d\t%d\t%d\t", chr, bundleno, record[s].pos, record[reale].pos);
        printf("%.2f\n", avg);
        bundleno++;
        
        for (int i = s; i < e; i++) {
            if (record[i].pos) {
                int pos = 0;
                if (strcmp(record[i].type, "jend") == 0) pos = 1;
                else if (strcmp(record[i].type, "jstart") == 0) pos = 2;
                double cov_val = get_next_cov(i, e);
                
                printf("%s\t%d\t", record[i].type, record[i].pos);
                printf("%.6f\t", record[i].change_perc);
                printf("%.0f\t", record[i].pos_cov);
                printf("%.3f", cov_val);
                
                if (pos) {
                    int nj_rec = record[i].num_indices;
                    for (int j = 0; j < nj_rec; j++) {
                        if (junc[record[i].indices[j]].cov > 0) {
                            int junc_pos = (pos == 1) ? junc[record[i].indices[j]].start : junc[record[i].indices[j]].end;
                            printf("\t%d:%c:%.0f", junc_pos, junc[record[i].indices[j]].strand, junc[record[i].indices[j]].cov);
                        }
                    }
                } else {
                    printf("\t%.2f", drop_arr[record[i].indices[0]].covdiff);
                }
                printf("\n");
            }
        }
    }
    
    return bundleno;
}

static int process_records(const char *chr, int bundleno) {
    int n = record_size;
    
    if (strcmp(record[n - 1].type, "tend") != 0) {
        exit(1);
    }
    
    int i = 1;
    int lasts = 0;
    int laste = 0;
    int lastjs = 0;
    int lastje = 0;
    int bundlend = record[0].pos;
    int s = 0;
    
    while (i < n) {
        if (strcmp(record[i].type, "tstart") == 0) {
            if (record[i].change_perc > 0) {
                if (lastje && record[i].pos - record[lastje].pos < smallwin && record[lastje].change_perc < 0.5) {
                    int j = i - 1;
                    while (j >= 0 && (record[j].pos == 0 || record[j].pos == record[i].pos)) {
                        if (record[i].pos == record[j].pos) {
                            record[i].pos_cov = 0;
                            break;
                        }
                        j--;
                    }
                    record[i].pos = 0;
                } else {
                    lasts = i;
                }
            } else {
                if (record[i - 1].pos == record[i].pos && strcmp(record[i - 1].type, "jend") == 0) {
                    record[i].pos = 0;
                    record[i].pos_cov = 0;
                } else if (record[i].pos > bundlend) {
                    bundleno = print_small_bundle(chr, bundleno, s, i);
                    s = i;
                }
            }
        } else if (strcmp(record[i].type, "tend") == 0) {
            if (record[i].change_perc > 0) {
                if (lastjs && record[i].pos - record[lastjs].pos < smallwin && record[lastjs].change_perc < 0.5) {
                    int j = i - 1;
                    while (j >= 0 && (record[j].pos == 0 || record[j].pos == record[i].pos)) {
                        if (record[i].pos == record[j].pos) {
                            record[i].pos_cov = 0;
                            break;
                        }
                        j--;
                    }
                    record[i].pos = 0;
                } else {
                    laste = i;
                }
            } else if (record[i - 1].pos == record[i].pos && strcmp(record[i - 1].type, "jstart") == 0) {
                record[i].pos = 0;
                record[i].pos_cov = 0;
            }
        } else if (strcmp(record[i].type, "jstart") == 0) {
            int nj_rec = record[i].num_indices;
            int found_valid = 0;
            for (int j = 0; j < nj_rec; j++) {
                if (junc[record[i].indices[j]].cov > 0) {
                    found_valid = 1;
                    if (junc[record[i].indices[j]].end > bundlend) {
                        bundlend = junc[record[i].indices[j]].end;
                    }
                }
            }
            if (found_valid) {
                if (record[i].change_perc < 0.5 && laste && record[i].pos - record[laste].pos < smallwin) {
                    int j = laste - 1;
                    while (j >= 0 && (record[j].pos == 0 || record[j].pos == record[laste].pos)) {
                        if (record[laste].pos == record[j].pos) {
                            record[laste].pos_cov = 0;
                            break;
                        }
                        j--;
                    }
                    record[laste].pos = 0;
                    laste = 0;
                }
                lastjs = i;
            } else {
                int j = i - 1;
                while (j >= 0 && (record[j].pos == 0 || record[j].pos == record[i].pos)) {
                    if (record[i].pos == record[j].pos) {
                        record[i].pos_cov = 0;
                        break;
                    }
                    j--;
                }
                record[i].pos = 0;
            }
        } else {
            int nj_rec = record[i].num_indices;
            int found_valid = 0;
            for (int j = 0; j < nj_rec; j++) {
                if (junc[record[i].indices[j]].cov > 0) {
                    if (record[i].change_perc < 0.5 && lasts && record[i].pos - record[lasts].pos < smallwin) {
                        int k = lasts - 1;
                        while (k >= 0 && (record[k].pos == 0 || record[k].pos == record[lasts].pos)) {
                            if (record[lasts].pos == record[k].pos) {
                                record[lasts].pos_cov = 0;
                                break;
                            }
                            k--;
                        }
                        record[lasts].pos = 0;
                        lasts = 0;
                    }
                    lastje = i;
                    found_valid = 1;
                    break;
                }
            }
            if (!found_valid) {
                record[i].pos = 0;
            }
        }
        i++;
    }
    
    bundleno = print_small_bundle(chr, bundleno, s, i);
    
    return bundleno;
}

static int process_bundle(const char *chr, int bundleno) {
    int nb = covg_size;
    if (nb == 0) return bundleno;
    
    double avg = 0;
    int len = 0;
    
    typedef struct {
        int si;
        int ei;
    } ContRegion;
    
    ContRegion *s_arr = NULL;
    int s_size = 0;
    int s_cap = 0;
    
    int preve = -1;
    int seengoodavg = 0;
    double runavg = 0;
    
    for (int i = 0; i < nb; i++) {
        if (covg[i].start - 1 > preve) {
            if (!seengoodavg && runavg > 0 && s_size > 0) {
                double region_len = covg[s_arr[s_size - 1].ei].end - covg[s_arr[s_size - 1].si].start + 1;
                if (runavg / region_len > lowcov) seengoodavg = 1;
            }
            runavg = 0;
            if (s_size >= s_cap) {
                s_cap = s_cap ? s_cap * 2 : 64;
                s_arr = realloc(s_arr, s_cap * sizeof(ContRegion));
            }
            s_arr[s_size].si = i;
            s_arr[s_size].ei = i;
            s_size++;
        }
        
        s_arr[s_size - 1].ei = i;
        preve = covg[i].end;
        int clen = covg[i].end - covg[i].start + 1;
        avg += covg[i].cov * clen;
        if (!seengoodavg) runavg += covg[i].cov * clen;
        len += clen;
    }
    
    avg /= len;
    if (runavg > 0 && s_size > 0) {
        double region_len = covg[s_arr[s_size - 1].ei].end - covg[s_arr[s_size - 1].si].start + 1;
        if (runavg / region_len > lowcov) seengoodavg = 1;
    }
    
    if (len > win && seengoodavg) {
        int nj = junc_size;
        sort_jend(nj);
        
        process_junctions(nj);
        
        clear_drop();
        
        int js = 0;
        int je = 0;
        for (int i = 0; i < s_size; i++) {
            get_drop(s_arr[i].si, s_arr[i].ei, nb, &js, &je, js, je, nj);
        }
        
        int nd = drop_size - 1;
        
        // Clear records
        for (int i = 0; i < record_size; i++) {
            free(record[i].indices);
        }
        record_size = 0;
        
        int ib = 0;
        int id = 0;
        js = 0;
        je = 0;
        
        int idx0 = 0;
        push_record("tstart", covg[0].start, &idx0, 1, 0, covg[0].cov, 0);
        id++;
        
        int prevpos = covg[0].start;
        
        while (id < nd || js < nj || je < nj) {
            get_record(&ib, &id, &js, &je, &prevpos, nb, nd, nj, chr);
        }
        
        if (covg[nb - 1].end > prevpos + 1) {
            int tmpib = ib;
            double avgcov = get_cov(prevpos + 1, covg[nb - 1].end - 1, &tmpib, nb);
            record[record_size - 1].cov_to_next += avgcov;
        }
        
        int idx_nd = nd;
        push_record("tend", covg[nb - 1].end, &idx_nd, 1, 0, covg[nb - 1].cov, 0);
        
        bundleno = process_records(chr, bundleno);
    }
    
    free(s_arr);
    clear_covg();
    clear_junc();
    
    return bundleno;
}

static int add_procjunc_to_bundle(int *bundleend, const char *chr) {
    int nj = unprocjunc_size;
    
    if (nj == 0) return 1;
    
    int j = 0;
    while (j < nj) {
        if (strcmp(unprocjunc[j].chrname, chr) != 0 || unprocjunc[j].start > *bundleend) {
            if (j > 0) {
                for (int i = 0; i < j; i++) {
                    push_junc_main(unprocjunc[i].chrname, unprocjunc[i].start, unprocjunc[i].end,
                                   unprocjunc[i].cov, unprocjunc[i].strand, unprocjunc[i].ps);
                }
                // Shift remaining
                for (int i = j; i < nj; i++) {
                    unprocjunc[i - j] = unprocjunc[i];
                }
                unprocjunc_size -= j;
            }
            return 0;
        }
        if (unprocjunc[j].end > *bundleend) {
            *bundleend = unprocjunc[j].end;
        }
        j++;
    }
    
    if (j > 0) {
        for (int i = 0; i < j; i++) {
            push_junc_main(unprocjunc[i].chrname, unprocjunc[i].start, unprocjunc[i].end,
                           unprocjunc[i].cov, unprocjunc[i].strand, unprocjunc[i].ps);
        }
        unprocjunc_size = 0;
    }
    return 1;
}

static int add_junc_to_bundle(const char *chr, int bundleend, FILE *fJ) {
    char line[4096];
    
    while (fgets(line, sizeof(line), fJ)) {
        char chrname[256];
        int start, end;
        char name[256];
        double cov_val;
        char strand;
        char percs[256];
        
        // Parse line
        char strand_str[8];
        if (sscanf(line, "%255s\t%d\t%d\t%255s\t%lf\t%7s\t%255s",
                   chrname, &start, &end, name, &cov_val, strand_str, percs) != 7) {
            continue;
        }
        strand = strand_str[0];
        
        // Parse percentages
        double ps, po, pl, pr;
        if (sscanf(percs, "%lf-%lf-%lf-%lf", &ps, &po, &pl, &pr) != 4) {
            continue;
        }
        
        double d = (ps < po) ? ps : po;
        double p = (pl < pr) ? pl : pr;
        
        int last = 0;
        if (strcmp(chrname, chr) != 0 || start > bundleend) {
            last = 1;
        }
        
        if ((ps > splicenoise || po > splicenoise) && 
            ((cov_val >= lowcov && strand != '.') || cov_val >= smallcov)) {
            
            if ((cov_val > highcov && strand != '.') || 
                (p > splicenoise && (d > splicenoise || (d > highnoise && cov_val > smallcov && strand != '.')))) {
                
                end++;  // Adjust end
                if (!last) {
                    if (end > bundleend) {
                        bundleend = end;
                    }
                    push_junc_main(chrname, start, end, cov_val, strand, ps);
                } else {
                    push_unprocjunc(chrname, start, end, cov_val, strand, ps);
                }
            }
        }
        
        if (last) return bundleend;
    }
    
    return bundleend;
}

int main(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <coverage.bedgraph> <junctions.bed>\n", argv[0]);
        return 1;
    }
    
    const char *covfile = argv[1];
    const char *juncfile = argv[2];
    
    FILE *C = fopen(covfile, "r");
    if (!C) {
        fprintf(stderr, "Cannot open coverage file: %s\n", covfile);
        return 1;
    }
    
    FILE *fJ = fopen(juncfile, "r");
    if (!fJ) {
        fprintf(stderr, "Cannot open junction file: %s\n", juncfile);
        fclose(C);
        return 1;
    }
    
    // Skip first lines (track lines)
    char line[4096];
    fgets(line, sizeof(line), C);
    fgets(line, sizeof(line), fJ);
    
    int bundleno = 0;
    char chr[256] = "";
    int bundleend = 0;
    
    while (fgets(line, sizeof(line), C)) {
        char chrname[256];
        int start, end;
        double cov_val;
        
        if (sscanf(line, "%255s\t%d\t%d\t%lf", chrname, &start, &end, &cov_val) != 4) {
            continue;
        }
        start++;  // Adjust to 1-based
        
        if (start > bundleend + 1 || strcmp(chrname, chr) != 0) {
            bundleno = process_bundle(chr, bundleno);
            if (strcmp(chr, chrname) != 0) {
                strncpy(chr, chrname, 255);
                chr[255] = '\0';
                fprintf(stderr, "Finding %s TSS/TES candidates\n", chr);
            }
            bundleend = 0;
        }
        
        if (end > bundleend) bundleend = end;
        push_covg(start, end, cov_val);
        
        int toadd = add_procjunc_to_bundle(&bundleend, chr);
        
        if (toadd) {
            bundleend = add_junc_to_bundle(chr, bundleend, fJ);
        }
    }
    
    fclose(C);
    fclose(fJ);
    
    // Process last bundle
    bundleno = process_bundle(chr, bundleno);
    
    // Cleanup
    free(covg);
    free(junc);
    free(unprocjunc);
    free(drop_arr);
    for (int i = 0; i < record_size; i++) {
        free(record[i].indices);
    }
    free(record);
    free(jend);
    
    return 0;
}