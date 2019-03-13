#include "driver_state.h"
#include <cstring>
#include <cfloat>
#include <vector>
using namespace std;
driver_state::driver_state()
{
}

driver_state::~driver_state()
{
    delete [] image_color;
    delete [] image_depth;
}

enum clip_case {none, A_in, B_in, C_in, AB_in, BC_in, AC_in, all};
enum axis {Xax, Yax, Zax, err};

// This function should allocate and initialize the arrays that store color and
// depth.  This is not done during the constructor since the width and height
// are not known when this class is constructed.
void initialize_render(driver_state& state, int width, int height)
{
    state.image_width=width;
    state.image_height=height;
    unsigned max = width*height; 
    state.image_color = new pixel[max]; //Sets the array to the resolution of the screen.
    state.image_depth = new float[max];
    for(size_t i = 0; i < max; i++){
        state.image_color[i] = make_pixel(0, 0, 0); //initializes all colors to black,
        state.image_depth[i] = FLT_MAX;
    }
}

// Takes in the two coordinate arrays and returns an array containing
// the maximum X coordinate and the maximum Y coordinate
void getMaxXY(float x[3], float y[3], float *maxes, int w, int h){
  
    float maxX = x[0];
    float maxY = y[0];
    for(unsigned i = 1; i < 3; i++){
        if(x[i] > maxX)
            maxX = x[i];
        if(y[i] > maxY)
            maxY = y[i];
    }
    // 0 = max X; 1 = max Y
    maxes[0] = maxX;
    if(maxX >= w)
        maxes[0] = w-1;
    maxes[1] = maxY;
    if(maxY >= h)
        maxes[1] = h-1;

    
}
// Takes in the two coordinate arrays and returns an array containing
// the minum X coordinate and the minimum Y coordinate
void getMinXY(float x[3], float y[3], float* mins){
    float minX = x[0];
    float minY = y[0];
    for(unsigned i = 1; i < 3; i++){
        if(x[i] < minX)
            minX = x[i];
        if(y[i] < minY)
            minY = y[i];
    }
    // [0] = min X, [1] = min Y
    mins[0] = minX;
    if(minX < 0)
        mins[0] = 0;
    mins[1] = minY;
    if(minY < 0)
        mins[1] = 0;

}

// This function will be called to render the data that has been stored in this class.
// Valid values of type are:
//   render_type::triangle - Each group of three vertices corresponds to a triangle.
//   render_type::indexed -  Each group of three indices in index_data corresponds
//                           to a triangle.  These numbers are indices into vertex_data.
//   render_type::fan -      The vertices are to be interpreted as a triangle fan.
//   render_type::strip -    The vertices are to be interpreted as a triangle strip.
void render(driver_state& state, render_type type)
{
    data_geometry *tri; //Stores the data_geometry for each triangle
    int index; //controls which vertex is accessed
    float *ptr; // pointer to data in vertex array
    int *in; //Index ptr for indexed triangles
    data_vertex *v; //array that holds vertex data
    
    switch(type){
        case render_type::triangle:
            tri = new data_geometry[3];
            v = new data_vertex[3];
            index = 0;
            ptr = state.vertex_data;
            for(unsigned i = 0; i <(unsigned) state.num_vertices; i++){
                tri[index].data = ptr;
                v[index].data = ptr;
                ptr += state.floats_per_vertex; //increment ptr to point to beginning of next vertex data
                state.vertex_shader(v[index], tri[index], state.uniform_data);
                if(index == 2){ //every 3 vertices is rasterized as a triangle.
                    index = -1;
                    clip_triangle(state,(const data_geometry**) &tri);
                }
                index++;
            }
            delete []tri; //Ensure no memory leaks
            delete []v; //Ensure no memory leaks
            break;
        case render_type::indexed:
            tri = new data_geometry[3];
            v = new data_vertex[3];
            index = 0;
            in = state.index_data;
            ptr = state.vertex_data;
            for(unsigned i = 0; i < state.num_triangles*3; i++){
                tri[index].data = &state.vertex_data[in[i]*state.floats_per_vertex];
                v[index].data = &state.vertex_data[in[i]*state.floats_per_vertex];
                state.vertex_shader(v[index], tri[index], state.uniform_data);
                if(index == 2){
                    index = -1;
                    clip_triangle(state, (const data_geometry**) &tri);
                }
                index++;
            }
            delete []tri;
            delete []v;            
            break;
        case render_type::fan:
            tri = new data_geometry[3];
            v = new data_vertex[3];
            index = 0;
            ptr = state.vertex_data;
            tri[0].data = ptr;
            v[0].data = ptr;
            ptr += state.floats_per_vertex;
            state.vertex_shader(v[0], tri[0], state.uniform_data);
            for(unsigned i = 1; i < state.num_vertices; i++){
                tri[1].data = ptr;
                tri[2].data = ptr + state.floats_per_vertex;
                v[1].data = ptr;
                v[2].data = ptr + state.floats_per_vertex;
                ptr += state.floats_per_vertex;
                state.vertex_shader(v[1], tri[1], state.uniform_data);
                state.vertex_shader(v[2], tri[2], state.uniform_data);
                clip_triangle(state, (const data_geometry**) &tri);
            }
            delete []tri;
            delete []v;
            break;
        case render_type::strip: 
            tri = new data_geometry[3];
            ptr = state.vertex_data;
            v = new data_vertex[3];
            index = 0;
            //Construct and rasterize initial triangle
            tri[0].data = ptr; //A
            v[0].data = ptr;
            ptr += state.floats_per_vertex;
            state.vertex_shader(v[0], tri[0], state.uniform_data);
            
            tri[1].data = ptr; //B
            v[1].data = ptr;
            ptr += state.floats_per_vertex;
            state.vertex_shader(v[1], tri[1], state.uniform_data);
            
            tri[2].data = ptr; //C
            v[2].data = ptr; 
            ptr += state.floats_per_vertex;
            state.vertex_shader(v[2], tri[2], state.uniform_data);
            clip_triangle(state, (const data_geometry**) &tri);

            for(unsigned i = 1; i <(unsigned) state.num_vertices -2; i++){
                if(i % 2){
                    tri[0].data = tri[1].data;
                    v[0].data = v[1].data;
                    tri[1].data = ptr;
                    v[1].data = ptr;
                    state.vertex_shader(v[0], tri[0], state.uniform_data);
                    state.vertex_shader(v[1], tri[1], state.uniform_data);
                    state.vertex_shader(v[2], tri[2], state.uniform_data);
                }
                else{
                    tri[0].data = tri[2].data;
                    v[0].data = tri[2].data;
                    tri[2].data = ptr;
                    v[2].data = ptr; 
                    state.vertex_shader(v[0], tri[0], state.uniform_data);
                    state.vertex_shader(v[1], tri[1], state.uniform_data);
                    state.vertex_shader(v[2], tri[2], state.uniform_data);
                }
                ptr += state.floats_per_vertex;
                clip_triangle(state, (const data_geometry**) &tri);
            }
            delete []tri; //Ensure no memory leaks
            delete []v; //Ensure no memory leaks
            break;
        default: //Error states
            break;
    }
}

axis getAxis(int face){
    if(face == 0 || face == 1)
        return Xax;
    if(face == 2 || face ==3)
        return Yax;
    if(face == 4 || face == 5)
        return Zax;
    return err;//Indicates rror
}

int getSign(int face){
    if(face%2 == 0)
        return 1;
    
    return -1;
        
}

clip_case getCase(const data_geometry* in[3], int sign, axis ax){
    enum clip_case clip = none; 
    bool Ain = false, Bin = false, Cin = false;
    float w[3];
    float pointA = (*in)[vA].gl_Position[ax], pointB = (*in)[vB].gl_Position[ax], 
        pointC = (*in)[vC].gl_Position[ax];
    for(int i = 0; i < 3;i++){
        w[i] = (*in)[i].gl_Position[3];
    }
    
    //Determine which vertices are within range
    if(sign > 0){
        if(pointA < w[vA])
            Ain = true;
        if(pointB < w[vB])
            Bin = true;
        if(pointC < w[vC])
            Cin = true;
    }
    else{
        if(pointA > -w[vA])
            Ain = true;
        if(pointB > -w[vB])
            Bin = true;
        if(pointC > -w[vC])
            Cin = true;
    }

    //Output proper case number
    if(!Ain && !Bin && !Cin){
        clip = none;
    }
    else if(Ain && !Bin && !Cin){
        clip = A_in;
    }
    else if(!Ain && Bin && !Cin){
        clip = B_in;
    }
    else if(!Ain && !Bin && Cin){
        clip = C_in;
    }
    else if(Ain && Bin && !Cin){
        clip = AB_in;
    }
    else if(!Ain && Bin && Cin){
        clip = BC_in;
    }
    else if(Ain && !Bin && Cin){
        clip = AC_in;
    }
    else if(Ain && Bin && Cin){
        clip = all;
    }
    
    return clip;    
}


void dataFillOne(driver_state& state, data_geometry dg[3], const data_geometry* in[3], float weight_1,
        float weight_2, vec4 p1, vec4 p2, vec4 p3, unsigned index0, unsigned index1, unsigned index2){

    float noPerWgt1;
    float noPerWgt2;
    float noPer1;
    float noPer2;
    for(int i = 0; i < MAX_FLOATS_PER_VERTEX; i++){
        switch(state.interp_rules[i]){
            case interp_type::flat:
                dg[0].data[i] = (*in)[index0].data[i];
                dg[1].data[i] = (*in)[index0].data[i];
                dg[2].data[i] = (*in)[index0].data[i];
                break;
            case interp_type::smooth:
                dg[0].data[i] = (*in)[index0].data[i];
                dg[1].data[i] = weight_1*(*in)[index0].data[i] + (1 - weight_1)*(*in)[index1].data[i];
                dg[2].data[i] = weight_2*(*in)[index0].data[i] + (1 - weight_2)*(*in)[index2].data[i];
                break;
            case interp_type::noperspective:
                noPerWgt1 = 1.0f / (weight_1*p1[3] + (1-weight_1)*p2[3]);
                noPerWgt2 = 1.0f / (weight_2*p1[3] + (1-weight_2)*p3[3]);
                noPer1 = weight_1*p1[3]*noPerWgt1;
                noPer2 = weight_2*p1[3]*noPerWgt2;
                dg[0].data[i] = (*in)[index0].data[i];
                dg[1].data[i] = noPer1*(*in)[index0].data[i] + (1 - noPer1)*(*in)[index1].data[i];
                dg[2].data[i] = noPer2*(*in)[index0].data[i] + (1 - noPer2)*(*in)[index2].data[i];
                break;
            default:
                break;
            }
        }
}  
void dataFillTwo(driver_state& state, data_geometry dg1[3], data_geometry dg2[3], const data_geometry* in[3],
 float weight_1, float weight_2, vec4 p1, vec4 p2, vec4 p3, unsigned index0, unsigned index1, unsigned index2){

    float noPerWgt1;
    float noPerWgt2;
    float noPer1;
    float noPer2;
    for(int i = 0; i < MAX_FLOATS_PER_VERTEX; i++){
        switch(state.interp_rules[i]){
            case interp_type::flat:
                dg1[0].data[i] = (*in)[index1].data[i];
                dg1[1].data[i] = (*in)[index1].data[i];
                dg1[2].data[i] = (*in)[index1].data[i];

                dg2[0].data[i] = (*in)[index2].data[i];
                dg2[1].data[i] = (*in)[index2].data[i];
                dg2[2].data[i] = (*in)[index2].data[i];
                break;
            case interp_type::smooth:
                dg1[0].data[i] = (*in)[index1].data[i];
                dg1[1].data[i] = (*in)[index2].data[i];
                dg1[2].data[i] = weight_2*(*in)[index2].data[i] + (1 - weight_2)*(*in)[index0].data[i];

                dg2[0].data[i] = (*in)[index1].data[i];
                dg2[1].data[i] = weight_2*(*in)[index2].data[i] + (1-weight_2)*(*in)[index0].data[i];
                dg2[2].data[i] = weight_1*(*in)[index1].data[i] + (1 - weight_1)*(*in)[index0].data[i];
                break;
            case interp_type::noperspective:
                noPerWgt1 = 1.0f / (weight_1*p2[3] + (1-weight_1)*p1[3]);
                noPerWgt2 = 1.0f / (weight_2*p3[3] + (1-weight_2)*p1[3]);
                noPer1 = weight_1*p2[3]*noPerWgt1;
                noPer2 = weight_2*p3[3]*noPerWgt2;
                dg1[0].data[i] = (*in)[index1].data[i];
                dg1[1].data[i] = (*in)[index2].data[i];
                dg1[2].data[i] = noPer2*(*in)[index2].data[i] + (1 - noPer2)*(*in)[index0].data[i];

                dg2[0].data[i] = (*in)[index1].data[i];
                dg2[1].data[i] = noPer2*(*in)[index2].data[i] + (1 - noPer2)*(*in)[index0].data[i];
                dg2[2].data[i] = noPer1*(*in)[index1].data[i] + (1 - noPer1)*(*in)[index0].data[i];
                break;
            default:
                break;
            }
        }
}  
void createTriangle(driver_state& state, vector<data_geometry*>& tri, const data_geometry* in[3],
    int sign, axis ax, unsigned index0, unsigned index1, unsigned index2, bool oneIn){
    vec4 p1 = (*in)[index0].gl_Position; 
    vec4 p2 = (*in)[index1].gl_Position; 
    vec4 p3 = (*in)[index2].gl_Position;
    float weight_1;
    float weight_2;
    vec4 newP2;
    vec4 newP3; 
    //Determine line weights to the outside points
    //Create new points at cut off plane
    //Create new triangle
    if(oneIn){
        weight_1 = (sign*p2[3] - p2[ax]) / (p1[ax] - sign*p1[3] + sign*p2[3] - p2[ax]); 
        weight_2 = (sign*p3[3] - p3[ax]) / (p1[ax] - sign*p1[3] + sign*p3[3] - p3[ax]);
        newP2 = weight_1 * p1 + (1 - weight_1)*p2;
        newP3 = weight_2 * p1 + (1 - weight_2)*p3;
        data_geometry* dg1 = new data_geometry[3];
        dg1[0].gl_Position = p1;
        dg1[0].data = new float[MAX_FLOATS_PER_VERTEX];
        dg1[1].gl_Position = newP2;
        dg1[1].data = new float[MAX_FLOATS_PER_VERTEX];
        dg1[2].gl_Position  = newP3;
        dg1[2].data = new float[MAX_FLOATS_PER_VERTEX];
        dataFillOne(state, dg1, in, weight_1, weight_2, p1, p2, p3, index0, index1, index2); 
        tri.push_back(dg1);     
    }
    else{   
        weight_1 = (sign*p1[3] - p1[ax]) / (p2[ax] - sign*p2[3] + sign*p1[3] - p1[ax]); 
        weight_2 = (sign*p1[3] - p1[ax]) / (p3[ax] - sign*p3[3] + sign*p1[3] - p1[ax]);
        newP2 = weight_1 * p2 + (1 - weight_1)*p1;
        newP3 = weight_2 * p3 + (1 - weight_2)*p1;
        data_geometry* dg1 = new data_geometry[3];
        data_geometry* dg2 = new data_geometry[3];
        //first triangle
        dg1[0].gl_Position = p2;
        dg1[0].data = new float[MAX_FLOATS_PER_VERTEX];
        dg1[1].gl_Position = p3;
        dg1[1].data = new float[MAX_FLOATS_PER_VERTEX];
        dg1[2].gl_Position  = newP3;
        dg1[2].data = new float[MAX_FLOATS_PER_VERTEX];
        //Second triangle
        dg2[0].gl_Position = p2;
        dg2[0].data = new float[MAX_FLOATS_PER_VERTEX];
        dg2[1].gl_Position = newP3;
        dg2[1].data = new float[MAX_FLOATS_PER_VERTEX];
        dg2[2].gl_Position  = newP2;
        dg2[2].data = new float[MAX_FLOATS_PER_VERTEX];
        dataFillTwo(state, dg1, dg2, in, weight_1, weight_2, p1, p2, p3, index0, index1, index2); 
        tri.push_back(dg1);     
        tri.push_back(dg2);
    }
    
}

void deleteTriData(data_geometry* in[3]){
    for(unsigned i = 0; i < 3; i++){
        delete[] (*in)[i].data;
    }

    delete in;
}
// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.
void clip_triangle(driver_state& state, const data_geometry* in[3],int face)
{
    vector<data_geometry*> tri;
    if(face==6)
    {
        rasterize_triangle(state, in);
        return;
    }
    else{
        int sign = getSign(face);
        enum axis ax = getAxis(face);
        enum clip_case clip = getCase(in, sign, ax);
        unsigned i;
        switch(clip){
            case none:
                return;
                break;
            case A_in:
                createTriangle(state, tri, in, sign, ax, vA, vB, vC, true);
                face++;
                for(i = 0; i < tri.size(); i++){
                    clip_triangle(state, (const data_geometry**) &tri.at(i), face); 
                }
                break;
            case B_in:
                createTriangle(state, tri, in, sign, ax, vB, vC, vA, true);
                face++;
                for(i = 0; i < tri.size(); i++){
                    clip_triangle(state, (const data_geometry**) &tri.at(i), face); 
                }
                break;
            case C_in:
                createTriangle(state, tri, in, sign, ax, vC, vA, vB, true);
                face++;
                for(i = 0; i < tri.size(); i++){
                    clip_triangle(state, (const data_geometry**) &tri.at(i), face); 
                }
                break;
            case AB_in:
                createTriangle(state, tri, in, sign, ax, vC, vA, vB, false);
                face++;
                for(i = 0; i < tri.size(); i++){
                    clip_triangle(state, (const data_geometry**) &tri.at(i), face); 
                }
                break;
            case BC_in:
                createTriangle(state, tri, in, sign, ax, vA, vB, vC, false);
                face++;
                for(i = 0; i < tri.size(); i++){
                    clip_triangle(state, (const data_geometry**) &tri.at(i), face); 
                }
                break;
            case AC_in:
                createTriangle(state, tri, in, sign, ax, vB, vC, vA, false);
                face++;
                for(i = 0; i < tri.size(); i++){
                    clip_triangle(state, (const data_geometry**) &tri.at(i), face);
                }
                break;
            case all:
                face++;
                clip_triangle(state, in, face);
                break;
            default:
                break;
        }
    }
  //  for(unsigned j = 0; j < 2; j++){
    //    deleteTriData( (data_geometry**) tri.at(j));
   // }
    
}

void dataFragmentFill(driver_state& state, data_fragment& in, const data_geometry* dg[3], 
        float alpha_p, float beta_p, float gamma_p){
    float wa = (*dg)[vA].gl_Position[3], wb =(*dg)[vB].gl_Position[3], wc = (*dg)[vC].gl_Position[3];
    float k = (alpha_p / wa) + (beta_p / wb) + (gamma_p / wc);
    float alpha = alpha_p / (wa*k), beta = beta_p / (wb*k), gamma = gamma_p / (wc*k);
    for(unsigned i =0; i < MAX_FLOATS_PER_VERTEX; i++){
        switch(state.interp_rules[i]){
            case interp_type::flat:
                in.data[i] = (*dg)[vA].data[i];
                break;
            case interp_type::smooth: 
                in.data[i] = alpha*(*dg)[vA].data[i] + beta*(*dg)[vB].data[i] + gamma*(*dg)[vC].data[i];
                break;
            case interp_type::noperspective:
                in.data[i] = alpha_p*(*dg)[vA].data[i] + beta_p*(*dg)[vB].data[i] + gamma_p*(*dg)[vC].data[i]; 
                break;
            default:
                break;
        }
    }
}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry* in[3])
{
    
    int w = state.image_width;
    int h = state.image_height; 
    float x_coord[3], y_coord[3], z_coord[3], w_coord[3];// Arrays that hold the NDC for the vertices
    float wOvr2 = 0.5f*w; // (w/2)
    float hOvr2 = 0.5f*h; // (h/2)
    float wNDC = wOvr2 - 0.5f; // (w/2 -(1/2))
    float hNDC = hOvr2 - 0.5f; // (h/2 - (1/2))
    float alpha_p = 0.0f, beta_p = 0.0f, gamma_p = 0.0f, areaABC = 0.0f;
    float k0 = 0.0f, k1 = 0.0f, k2 = 0.0f; //For calculated the area of alpha
    float t0 = 0.0f, t1 = 0.0f, t2 = 0.0f; //For calculating the area of beta
    float z0 = 0.0f, z1 = 0.0f, z2 = 0.0f; //For calclating the area of gamma
    float mins[2], maxes[2];
    //Part 3, section 3
    for(unsigned k = 0; k < 3; k++){
        x_coord[k] = wOvr2*((*in)[k].gl_Position[0] / (*in)[k].gl_Position[3]) + wNDC;
        y_coord[k] = hOvr2*((*in)[k].gl_Position[1] / (*in)[k].gl_Position[3]) + hNDC;
        z_coord[k] = hOvr2*((*in)[k].gl_Position[2] / (*in)[k].gl_Position[3]) + hNDC;
        w_coord[k] = hOvr2*((*in)[k].gl_Position[3] / (*in)[k].gl_Position[3]) + hNDC;
    } 

    // vA = vertex A; vB = vertex B; vC = vertex C;
    //Part3, section 4
    areaABC = (x_coord[vB]*y_coord[vC] - x_coord[vC]*y_coord[vB]);
    areaABC -= (x_coord[vA]*y_coord[vC] - x_coord[vC]*y_coord[vA]);
    areaABC += (x_coord[vA]*y_coord[vB] - x_coord[vB]*y_coord[vA]);
    areaABC *= 0.5f;
    //alpha area
    k0 = (x_coord[vB]*y_coord[vC] - x_coord[vC]*y_coord[vB]);
    k1 = y_coord[vB] - y_coord[vC];
    k2 = x_coord[vC] - x_coord[vB];
    //beta area
    t0 = (x_coord[vC]*y_coord[vA] - x_coord[vA]*y_coord[vC]);
    t1 = y_coord[vC] - y_coord[vA];
    t2 = x_coord[vA] - x_coord[vC];
    //gamma area
    z0 = (x_coord[vA]*y_coord[vB] - x_coord[vB]*y_coord[vA]);
    z1 = y_coord[vA] - y_coord[vB];
    z2 = x_coord[vB] - x_coord[vA];
    getMinXY(x_coord, y_coord, mins);
    getMaxXY(x_coord, y_coord, maxes, state.image_width, state.image_height);
    for(int j = mins[1]+1; j < maxes[1] +1; j++){
        for(int i = mins[0]+1; i < maxes[0] + 1; i++){
            alpha_p = (0.5f*(k0 + k1*i + k2*j))/areaABC;
            beta_p = (0.5f*(t0 + t1*i + t2*j))/areaABC;
            gamma_p = (0.5f*(z0 + z1*i + z2*j))/areaABC;
            if(alpha_p >= 0 && beta_p >= 0 && gamma_p >=0){
                float z_val = alpha_p*z_coord[vA] + beta_p*z_coord[vB] + gamma_p*z_coord[vC];
                if(z_val < state.image_depth[i+j*w]){  
                    float *fragData = new float[MAX_FLOATS_PER_VERTEX];
                    state.image_depth[i+j*w] = z_val;
                    data_fragment df{fragData};
                    data_output out;
                    dataFragmentFill(state, df, in, alpha_p, beta_p, gamma_p);
                    state.fragment_shader(df, out, state.uniform_data);
                    state.image_color[i+j*w] =
                        make_pixel(out.output_color[0]*255,out.output_color[1]*255.0, out.output_color[2]*255.0);
                    delete[]  fragData;
                }
            }
        }
    }
}

