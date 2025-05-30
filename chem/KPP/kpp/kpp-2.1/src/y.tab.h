#define JACOBIAN 257
#define DOUBLE 258
#define FUNCTION 259
#define DEFVAR 260
#define DEFRAD 261
#define DEFFIX 262
#define SETVAR 263
#define SETRAD 264
#define SETFIX 265
#define HESSIAN 266
#define STOICMAT 267
#define STOCHASTIC 268
#define INITVALUES 269
#define EQUATIONS 270
#define LUMP 271
#define INIEQUAL 272
#define EQNEQUAL 273
#define EQNCOLON 274
#define LMPCOLON 275
#define LMPPLUS 276
#define SPCPLUS 277
#define SPCEQUAL 278
#define ATOMDECL 279
#define CHECK 280
#define CHECKALL 281
#define REORDER 282
#define MEX 283
#define DUMMYINDEX 284
#define EQNTAGS 285
#define LOOKAT 286
#define LOOKATALL 287
#define TRANSPORT 288
#define TRANSPORTALL 289
#define MONITOR 290
#define USES 291
#define SPARSEDATA 292
#define WRFCONFORM 293
#define WRITE_ATM 294
#define WRITE_SPC 295
#define WRITE_MAT 296
#define WRITE_OPT 297
#define INITIALIZE 298
#define XGRID 299
#define YGRID 300
#define ZGRID 301
#define USE 302
#define LANGUAGE 303
#define INTFILE 304
#define DRIVER 305
#define RUN 306
#define INLINE 307
#define ENDINLINE 308
#define PARAMETER 309
#define SPCSPC 310
#define INISPC 311
#define INIVALUE 312
#define EQNSPC 313
#define EQNSIGN 314
#define EQNCOEF 315
#define RATE 316
#define LMPSPC 317
#define SPCNR 318
#define ATOMID 319
#define LKTID 320
#define MNIID 321
#define INLCTX 322
#define INCODE 323
#define SSPID 324
#define EQNLESS 325
#define EQNTAG 326
#define EQNGREATER 327
#define TPTID 328
#define USEID 329
#ifdef YYSTYPE
#undef  YYSTYPE_IS_DECLARED
#define YYSTYPE_IS_DECLARED 1
#endif
#ifndef YYSTYPE_IS_DECLARED
#define YYSTYPE_IS_DECLARED 1
typedef union{
  char str[500];
} YYSTYPE;
#endif /* !YYSTYPE_IS_DECLARED */
extern YYSTYPE yylval;
