/*
 * These are the "generic" xdr routines.
 */
extern bool_t xdr_void (void);
extern bool_t xdr_short (XDR *__xdrs, short *__sp);
extern bool_t xdr_u_short (XDR *__xdrs, u_short *__usp);
extern bool_t xdr_int (XDR *__xdrs, int *__ip);
extern bool_t xdr_u_int (XDR *__xdrs, u_int *__up);
extern bool_t xdr_long (XDR *__xdrs, long *__lp);
extern bool_t xdr_u_long (XDR *__xdrs, u_long *__ulp);
extern bool_t xdr_hyper (XDR *__xdrs, quad_t *__llp);
extern bool_t xdr_u_hyper (XDR *__xdrs, u_quad_t *__ullp);
extern bool_t xdr_longlong_t (XDR *__xdrs, quad_t *__llp);
extern bool_t xdr_u_longlong_t (XDR *__xdrs, u_quad_t *__ullp);
extern bool_t xdr_int8_t (XDR *__xdrs, int8_t *__ip);
extern bool_t xdr_uint8_t (XDR *__xdrs, uint8_t *__up);
extern bool_t xdr_int16_t (XDR *__xdrs, int16_t *__ip);
extern bool_t xdr_uint16_t (XDR *__xdrs, uint16_t *__up);
extern bool_t xdr_int32_t (XDR *__xdrs, int32_t *__ip);
extern bool_t xdr_uint32_t (XDR *__xdrs, uint32_t *__up);
extern bool_t xdr_int64_t (XDR *__xdrs, int64_t *__ip);
extern bool_t xdr_uint64_t (XDR *__xdrs, uint64_t *__up);
extern bool_t xdr_quad_t (XDR *__xdrs, quad_t *__ip);
extern bool_t xdr_u_quad_t (XDR *__xdrs, u_quad_t *__up);
extern bool_t xdr_bool (XDR *__xdrs, bool_t *__bp);
extern bool_t xdr_enum (XDR *__xdrs, enum_t *__ep);
extern bool_t xdr_array (XDR * _xdrs, caddr_t *__addrp, u_int *__sizep,u_int __maxsize, u_int __elsize, xdrproc_t __elproc);
extern bool_t xdr_bytes (XDR *__xdrs, char **__cpp, u_int *__sizep,u_int __maxsize);
extern bool_t xdr_opaque (XDR *__xdrs, caddr_t __cp, u_int __cnt);
extern bool_t xdr_string (XDR *__xdrs, char **__cpp, u_int __maxsize);
extern bool_t xdr_union (XDR *__xdrs, enum_t *__dscmp, char *__unp, __const struct xdr_discrim *__choices,xdrproc_t __dfault);
extern bool_t xdr_char (XDR *__xdrs, char *__cp);
extern bool_t xdr_u_char (XDR *__xdrs, u_char *__cp);
extern bool_t xdr_vector (XDR *__xdrs, char *__basep, u_int __nelem,u_int __elemsize, xdrproc_t __xdr_elem);
extern bool_t xdr_float (XDR *__xdrs, float *__fp);
extern bool_t xdr_double (XDR *__xdrs, double *__dp);
extern bool_t xdr_reference (XDR *__xdrs, caddr_t *__xpp, u_int __size,xdrproc_t __proc);
extern bool_t xdr_pointer (XDR *__xdrs, char **__objpp,u_int __obj_size, xdrproc_t __xdr_obj);
extern bool_t xdr_wrapstring (XDR *__xdrs, char **__cpp);
extern u_long xdr_sizeof (xdrproc_t, void *);

extern bool_t xdrrec_endofrecord (XDR *__xdrs, bool_t __sendnow);
extern bool_t xdrrec_skiprecord (XDR *__xdrs);