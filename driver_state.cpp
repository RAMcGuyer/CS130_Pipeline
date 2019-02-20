#include "driver_state.h"
#include <cstring>

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
    state.image_color= new pixel[max]; //Sets the array to the resolution of the screen.
    for(size_t i = 0; i < max; i++){
        state.image_color[i] = make_pixel(0, 0, 0); //initializes all colors to black,
    }
    state.image_depth=0;
    std::cout<<"initialize state.image_depth."<<std::endl;
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
    data_geometry *tri;
    int index;
    float *ptr;
    data_vertex *v;
    switch(type){
        case render_type::triangle:
            tri = new data_geometry[3];
            v = new data_vertex[3];
            index = 0;
            ptr = state.vertex_data;
            for(unsigned i = 0; i <(unsigned) state.num_vertices; i++){
                tri[index].data = ptr;
                v[index].data = ptr;
                ptr += state.floats_per_vertex;
                if(index == 2){
                    index = -1;
                    state.vertex_shader(v[0], tri[0], state.uniform_data);
                    state.vertex_shader(v[1], tri[1], state.uniform_data);
                    state.vertex_shader(v[2], tri[2], state.uniform_data);
                    rasterize_triangle(state,(const data_geometry**) &tri);
                }
                index++;
            }
            delete []tri;
            delete []v;
            break;
        case render_type::indexed:
            
            break;
        case render_type::fan:
            
            break;
        case render_type::strip:
            
            break;
        default:
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

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry* in[3])
{
    
    auto w = state.image_width;
    auto h = state.image_height; 
    //Part 3, section 3

    
    int aX = 0, bX = 0, cX = 0;
    int aY = 0, bY = 0, cY = 0;
    aX = (w/2)*(*in)[0].gl_Position[0] + (w/2 -(1/2));
    bX = (w/2)*(*in)[1].gl_Position[0] + (w/2 -(1/2));
    cX = (w/2)*(*in)[2].gl_Position[0] + (w/2 -(1/2));
    
    aY = (h/2)*(*in)[0].gl_Position[1] + (h/2 -(1/2));
    bY = (h/2)*(*in)[1].gl_Position[1] + (h/2 -(1/2));
    cY = (h/2)*(*in)[2].gl_Position[1] + (h/2 -(1/2));

    //Part3, section 4
   
    //alpha area
    float k0 = (bX*cY - cX*bY);
    float k1 = bY - cY;
    float k2 = cX - bX;
    //beta area
    float t0 = (cX*aY-aX*cY);
    float t1 = cY - aY;
    float t2 = aX - cX;
    //gamma area
    float z0 = (aX*bY-bX*aY);
    float z1 = aY - bY;
    float z2 = bX - aX;

    float areaABC = 0.5f*((bX*cY - cX*bY) - (aX*cY - cX*aY) - (aX*bY - bX*aY));
    float alpha = 0, beta = 0, gamma =0;
    for(unsigned j = 0; j < h; j++){
        for(unsigned i = 0; i < w; i++){
            alpha = (0.5f*(k0 + k1*i + k2*j))/areaABC;
            beta = (0.5f*(t0 + t1*i + t2*j))/areaABC;
            gamma = (0.5f*(z0 + z1*i + z2*j))/areaABC;
            if(alpha >= 0 && beta >= 0 && gamma >=0){
                state.image_color[i+j*w] = make_pixel(255.0, 255.0, 255.0);
            }
        }
    }
}

