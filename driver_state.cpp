#include "driver_state.h"
#include <cstring>
#include <cfloat>

driver_state::driver_state()
{
}

driver_state::~driver_state()
{
    delete [] image_color;
    delete [] image_depth;
}

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
void getMaxXY(int x[3], int y[3], int *maxes, int w, int h){
  
    int maxX = x[0];
    int maxY = y[0];
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
void getMinXY(int x[3], int y[3], int* mins){
    int minX = x[0];
    int minY = y[0];
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
                    rasterize_triangle(state,(const data_geometry**) &tri);
                }
                index++;
            }
            delete []tri; //Ensure no memory leaks
            delete []v; //Ensure no memory leaks
            break;
        case render_type::indexed:
            
            break;
        case render_type::fan:
            
            break;
        case render_type::strip:
            
            break;
        default: //Error states
            break;
    }
    std::cout<<"TODO: implement rendering."<<std::endl;
}


// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.
void clip_triangle(driver_state& state, const data_geometry* in[3],int face)
{
    if(face==6)
    {
        rasterize_triangle(state, in);
        return;
    }
    std::cout<<"TODO: implement clipping. (The current code passes the triangle through without clipping them.)"<<std::endl;
    clip_triangle(state,in,face+1);
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
    
    auto w = state.image_width;
    auto h = state.image_height; 
    int x_coord[3], y_coord[3], z_coord[3], w_coord[3];// Arrays that hold the NDC for the vertices
    float wOvr2 = 0.5f*w; // (w/2)
    float hOvr2 = 0.5f*h; // (h/2)
    float wNDC = wOvr2 - 0.5f; // (w/2 -(1/2))
    float hNDC = hOvr2 - 0.5f; // (h/2 - (1/2))
    float alpha_p = 0.0f, beta_p = 0.0f, gamma_p = 0.0f, areaABC = 0.0f;
    float k0 = 0.0f, k1 = 0.0f, k2 = 0.0f; //For calculated the area of alpha
    float t0 = 0.0f, t1 = 0.0f, t2 = 0.0f; //For calculating the area of beta
    float z0 = 0.0f, z1 = 0.0f, z2 = 0.0f; //For calclating the area of gamma
    int mins[2], maxes[2];
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
    for(int j = mins[1]; j < maxes[1] +1; j++){
        for(int i = mins[0]; i < maxes[0] + 1; i++){
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

